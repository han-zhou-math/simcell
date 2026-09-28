#!/usr/bin/env python3
"""Visualize SimCELL snapshots written by ``outuvpx`` and the final solvers."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import SymLogNorm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("Data"))
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--nx", type=int, default=128)
    parser.add_argument("--ny", type=int, default=128)
    parser.add_argument("--nmarkers", type=int, default=100)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--output-every", type=int, default=1)
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--x-half", type=float, default=1.0)
    parser.add_argument("--y-half", type=float, default=1.0)
    return parser.parse_args()


def read_field(path: Path, shape: tuple[int, int]) -> np.ndarray:
    values = np.fromfile(path, dtype=np.float64)
    expected = int(np.prod(shape))
    if values.size != expected:
        raise ValueError(f"{path} contains {values.size} doubles; expected {expected}")
    return values.reshape(shape, order="F")


def read_frame(data_dir: Path, frame: int, nx: int, ny: int):
    u_face = read_field(data_dir / f"frun.u.{frame:04d}", (nx, ny))
    v_raw = read_field(data_dir / f"frun.v.{frame:04d}", (nx, ny - 1))
    pressure = read_field(data_dir / f"frun.p.{frame:04d}", (nx, ny))
    interface = np.loadtxt(data_dir / f"frun.ib.{frame:04d}")

    # Average the staggered face velocities to pressure-cell centers.  The
    # x direction is periodic; the y-wall values are zero.
    u_center = 0.5 * (u_face + np.roll(u_face, 1, axis=0))
    v_face = np.zeros((nx, ny + 1), dtype=np.float64)
    v_face[:, 1:ny] = v_raw
    v_center = 0.5 * (v_face[:, :-1] + v_face[:, 1:])
    return u_center, v_center, pressure, interface


def closed_curve(interface: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    closed = np.vstack([interface, interface[0]])
    return closed[:, 0], closed[:, 1]


def interface_metrics(interface: np.ndarray, h: float) -> dict[str, object]:
    segments = np.linalg.norm(np.roll(interface, -1, axis=0) - interface, axis=1)
    x, y = interface.T
    area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    return {
        "centroid": [float(x.mean()), float(y.mean())],
        "area": float(area),
        "spacing_min": float(segments.min()),
        "spacing_mean": float(segments.mean()),
        "spacing_max": float(segments.max()),
        "mean_spacing_over_bulk_h": float(segments.mean() / h),
    }


def parse_diagnostics(log_path: Path) -> dict[str, np.ndarray | float]:
    chemical = []
    actin = []
    jump = []
    brinkman = []
    initial_actin_mass = np.nan
    gmres_converged_count = 0
    geometry_skip_warning_count = 0
    mpi_bind_warning_count = 0
    for line in log_path.read_text(errors="replace").splitlines():
        gmres_converged_count += "GMRES converged." in line
        geometry_skip_warning_count += "special case:" in line
        mpi_bind_warning_count += "bind() failed:" in line
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "STAGE07_DUALCHEM_MASS":
            chemical.append((int(fields[1]), float(fields[2]), float(fields[3])))
        elif fields[0] == "STAGE07_PHYSICAL_JUMP":
            jump.append((int(fields[1]), float(fields[2]), float(fields[3]), float(fields[4])))
        elif fields[0] == "STAGE12_ACTIN_INITIAL_MASS":
            initial_actin_mass = float(fields[3])
        elif fields[0] == "STAGE12_ACTIN_MASS":
            actin.append((int(fields[1]), float(fields[2]), float(fields[3])))
        elif fields[0] == "STAGE13_BRINKMAN":
            brinkman.append(tuple(map(int, fields[1:5])))
    return {
        "chemical": np.asarray(chemical, dtype=float),
        "actin": np.asarray(actin, dtype=float),
        "jump": np.asarray(jump, dtype=float),
        "brinkman": np.asarray(brinkman, dtype=int),
        "initial_actin_mass": float(initial_actin_mass),
        "gmres_converged_count": gmres_converged_count,
        "geometry_skip_warning_count": geometry_skip_warning_count,
        "mpi_bind_warning_count": mpi_bind_warning_count,
    }


def set_spatial_axis(axis: plt.Axes, title: str, x_half: float, y_half: float) -> None:
    axis.set_title(title)
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_xlim(-x_half, x_half)
    axis.set_ylim(-y_half, y_half)
    axis.set_aspect("equal")


def add_colorbar(fig: plt.Figure, image, axis: plt.Axes, label: str) -> None:
    colorbar = fig.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    colorbar.set_label(label)


def add_velocity_streamlines(
    axis: plt.Axes,
    u: np.ndarray,
    v: np.ndarray,
    *,
    density: float = 1.2,
    x_half: float = 1.0,
    y_half: float = 1.0,
):
    """Overlay instantaneous lab-frame streamlines on a cell-centered field."""
    nx, ny = u.shape
    dx = 2.0 * x_half / nx
    dy = 2.0 * y_half / ny
    x = np.linspace(-x_half + 0.5 * dx, x_half - 0.5 * dx, nx)
    y = np.linspace(-y_half + 0.5 * dy, y_half - 0.5 * dy, ny)
    speed = np.hypot(u, v)
    speed_scale = max(float(speed.max()), np.finfo(float).tiny)
    linewidth = 0.45 + 1.15 * speed.T / speed_scale
    return axis.streamplot(
        x,
        y,
        u.T,
        v.T,
        density=density,
        color="white",
        linewidth=linewidth,
        arrowsize=0.7,
        integration_direction="both",
        maxlength=2.5,
        zorder=4,
    )


def make_movie(
    output: Path,
    frames: list[int],
    times: np.ndarray,
    u_all: np.ndarray,
    v_all: np.ndarray,
    p_all: np.ndarray,
    chemical_all: np.ndarray | None,
    interfaces: list[np.ndarray],
    fps: int,
    x_half: float,
    y_half: float,
) -> None:
    speed_all = np.hypot(u_all, v_all)
    speed_max = float(speed_all.max())
    pressure_centered = p_all - p_all.mean(axis=(1, 2), keepdims=True)
    pressure_max = float(np.max(np.abs(pressure_centered)))
    pressure_norm = SymLogNorm(
        linthresh=max(pressure_max * 1.0e-3, 1.0e-8),
        linscale=0.8,
        vmin=-pressure_max,
        vmax=pressure_max,
    )

    panel_count = 3 if chemical_all is not None else 2
    fig, axes = plt.subplots(
        1, panel_count, figsize=(6.4 * panel_count, 5.8), constrained_layout=True
    )
    fig.suptitle("SimCELL solution evolution")
    extent = (-x_half, x_half, -y_half, y_half)
    speed_image = axes[0].imshow(
        speed_all[0].T,
        origin="lower",
        extent=extent,
        cmap="magma",
        vmin=0.0,
        vmax=speed_max,
        interpolation="bilinear",
    )
    pressure_image = axes[1].imshow(
        pressure_centered[0].T,
        origin="lower",
        extent=extent,
        cmap="coolwarm",
        norm=pressure_norm,
        interpolation="bilinear",
    )
    set_spatial_axis(axes[0], "Velocity magnitude", x_half, y_half)
    set_spatial_axis(axes[1], "Pressure (spatial mean removed)", x_half, y_half)
    add_colorbar(fig, speed_image, axes[0], r"$|\mathbf{u}|$")
    add_colorbar(fig, pressure_image, axes[1], r"$p-\langle p\rangle$")

    chemical_image = None
    if chemical_all is not None:
        chemical_perturbation = chemical_all - 1.0
        chemical_limit = float(np.max(np.abs(chemical_perturbation)))
        chemical_image = axes[2].imshow(
            chemical_perturbation[0].T,
            origin="lower",
            extent=extent,
            cmap="coolwarm",
            vmin=-chemical_limit,
            vmax=chemical_limit,
            interpolation="bilinear",
        )
        set_spatial_axis(axes[2], "Chemical perturbation", x_half, y_half)
        add_colorbar(fig, chemical_image, axes[2], r"$c-1$")

    ix, iy = closed_curve(interfaces[0])
    interface_lines = [
        axes[0].plot(ix, iy, color="cyan", linewidth=2.0, zorder=5)[0],
        axes[1].plot(ix, iy, color="black", linewidth=2.0, zorder=5)[0],
    ]
    if chemical_all is not None:
        interface_lines.append(
            axes[2].plot(ix, iy, color="black", linewidth=2.0, zorder=5)[0]
        )
    patches_before_streamlines = set(axes[0].patches)
    streamlines = add_velocity_streamlines(
        axes[0], u_all[0], v_all[0], density=1.0,
        x_half=x_half, y_half=y_half
    )
    streamline_arrows = [
        patch for patch in axes[0].patches if patch not in patches_before_streamlines
    ]
    time_label = fig.text(0.5, 0.015, f"frame {frames[0]:d}   t = {times[0]:.2f}", ha="center")

    writer = FFMpegWriter(
        fps=fps,
        metadata={"title": f"SimCELL solution movie through t={times[-1]:.3g}"},
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart"],
    )
    with writer.saving(fig, output, dpi=140):
        for index, frame in enumerate(frames):
            speed_image.set_data(speed_all[index].T)
            pressure_image.set_data(pressure_centered[index].T)
            if chemical_image is not None:
                chemical_image.set_data((chemical_all[index] - 1.0).T)
            ix, iy = closed_curve(interfaces[index])
            for line in interface_lines:
                line.set_data(ix, iy)
            streamlines.lines.remove()
            for arrow in streamline_arrows:
                arrow.remove()
            patches_before_streamlines = set(axes[0].patches)
            streamlines = add_velocity_streamlines(
                axes[0], u_all[index], v_all[index], density=1.0,
                x_half=x_half, y_half=y_half
            )
            streamline_arrows = [
                patch
                for patch in axes[0].patches
                if patch not in patches_before_streamlines
            ]
            time_label.set_text(f"frame {frame:d}   t = {times[index]:.2f}")
            writer.grab_frame()
    plt.close(fig)


def make_final_figure(
    output: Path,
    u: np.ndarray,
    v: np.ndarray,
    pressure: np.ndarray,
    interface: np.ndarray,
    chemical: np.ndarray,
    network: np.ndarray,
    free: np.ndarray,
    jump: np.ndarray,
    final_time: float,
    x_half: float,
    y_half: float,
    case_title: str,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.8), constrained_layout=True)
    fig.suptitle(rf"{case_title} at $t={final_time:.2f}$")
    extent = (-x_half, x_half, -y_half, y_half)
    ix, iy = closed_curve(interface)
    spatial_fields = [
        (
            np.hypot(u, v),
            "magma",
            r"Velocity field: speed and streamlines",
            r"$|\mathbf{u}|$",
        ),
        (pressure - pressure.mean(), "coolwarm", "Pressure (mean removed)", r"$p-\langle p\rangle$"),
        (chemical - 1.0, "coolwarm", "Chemical perturbation", r"$c-1$"),
        (np.ma.masked_where(network == 0.0, network), "viridis", "Network actin", r"$\theta_n$"),
        (np.ma.masked_where(free == 0.0, free), "cividis", "Free actin", r"$\theta_c$"),
    ]
    for axis, (field, cmap, title, label) in zip(axes.flat[:5], spatial_fields):
        if title in {"Pressure (mean removed)", "Chemical perturbation"}:
            limit = float(np.max(np.abs(field)))
            image = axis.imshow(
                field.T,
                origin="lower",
                extent=extent,
                cmap=cmap,
                vmin=-limit,
                vmax=limit,
                interpolation="bilinear",
            )
        else:
            image = axis.imshow(field.T, origin="lower", extent=extent, cmap=cmap, interpolation="bilinear")
        axis.plot(
            ix,
            iy,
            color="white" if title != "Chemical perturbation" else "black",
            linewidth=1.5,
            zorder=5,
        )
        set_spatial_axis(axis, title, x_half, y_half)
        add_colorbar(fig, image, axis, label)
        if title == "Velocity field: speed and streamlines":
            add_velocity_streamlines(
                axis, u, v, density=1.25, x_half=x_half, y_half=y_half
            )

    marker_fraction = np.arange(jump.size) / jump.size
    axes[1, 2].plot(marker_fraction, jump, color="tab:purple", linewidth=2.0)
    axes[1, 2].axhline(0.0, color="0.3", linewidth=0.8)
    axes[1, 2].set_title("Physical concentration jump")
    axes[1, 2].set_xlabel("normalized interface coordinate")
    axes[1, 2].set_ylabel(r"$c_i-c_e$")
    axes[1, 2].grid(alpha=0.2)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def make_diagnostics_figure(output: Path, diagnostics: dict[str, object]) -> None:
    chemical = diagnostics["chemical"]
    actin = diagnostics["actin"]
    jump = diagnostics["jump"]
    initial_actin_mass = float(diagnostics["initial_actin_mass"])

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.0), constrained_layout=True)
    chemical_reference_mass = float(chemical[0, 2])
    axes[0].plot(chemical[:, 1], chemical[:, 2], linewidth=2.0)
    axes[0].axhline(
        chemical_reference_mass, color="0.4", linestyle="--", linewidth=1.0
    )
    axes[0].set_title("Chemical mass")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel("mass")

    axes[1].plot(actin[:, 1], actin[:, 2], linewidth=2.0)
    if np.isfinite(initial_actin_mass):
        axes[1].scatter([0.0], [initial_actin_mass], color="black", s=20, zorder=3)
    axes[1].set_title("Total actin mass")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("mass")

    axes[2].fill_between(jump[:, 1], jump[:, 2], jump[:, 3], alpha=0.25)
    axes[2].plot(jump[:, 1], jump[:, 2], linewidth=1.5, label="minimum")
    axes[2].plot(jump[:, 1], jump[:, 3], linewidth=1.5, label="maximum")
    axes[2].axhline(0.0, color="0.4", linewidth=0.8)
    axes[2].set_title("Interface concentration jump range")
    axes[2].set_xlabel("time")
    axes[2].set_ylabel(r"$c_i-c_e$")
    axes[2].legend(frameon=False)
    for axis in axes:
        axis.grid(alpha=0.2)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pattern = re.compile(r"frun\.ib\.(\d{4})$")
    frames = sorted(
        int(match.group(1))
        for path in args.data_dir.glob("frun.ib.*")
        if (match := pattern.search(path.name))
    )
    if not frames:
        raise RuntimeError(f"no interface snapshots found under {args.data_dir}")

    data = [read_frame(args.data_dir, frame, args.nx, args.ny) for frame in frames]
    u_all = np.stack([item[0] for item in data])
    v_all = np.stack([item[1] for item in data])
    p_all = np.stack([item[2] for item in data])
    interfaces = [item[3] for item in data]
    arrays = [u_all, v_all, p_all, *interfaces]
    if not all(np.isfinite(array).all() for array in arrays):
        raise ValueError("non-finite value found in time-dependent output")
    if any(interface.shape != (args.nmarkers, 2) for interface in interfaces):
        raise ValueError("unexpected interface snapshot shape")
    chemical_frame_paths = [args.data_dir / f"frun.c.{frame:04d}" for frame in frames]
    chemical_all = None
    if all(path.exists() for path in chemical_frame_paths):
        chemical_all = np.stack(
            [read_field(path, (args.nx, args.ny)) for path in chemical_frame_paths]
        )
        if not np.isfinite(chemical_all).all():
            raise ValueError("non-finite value found in chemical time series")

    chemical = read_field(args.data_dir / "stage07.chemical.final.bin", (args.nx, args.ny))
    network = read_field(args.data_dir / "stage12.network.final.bin", (args.nx, args.ny))
    free = read_field(args.data_dir / "stage12.free.final.bin", (args.nx, args.ny))
    jump = np.fromfile(args.data_dir / "stage07.jump.final.bin", dtype=np.float64)
    final_arrays = [chemical, network, free, jump]
    if not all(np.isfinite(array).all() for array in final_arrays):
        raise ValueError("non-finite value found in final coupled fields")

    times = np.asarray(frames, dtype=float) * args.dt * args.output_every
    diagnostics = parse_diagnostics(args.log)
    movie_path = args.output_dir / "solution_movie.mp4"
    final_path = args.output_dir / "solution_final.png"
    diagnostics_path = args.output_dir / "diagnostics.png"
    make_movie(
        movie_path,
        frames,
        times,
        u_all,
        v_all,
        p_all,
        chemical_all,
        interfaces,
        args.fps,
        args.x_half,
        args.y_half,
    )
    make_final_figure(
        final_path,
        u_all[-1],
        v_all[-1],
        p_all[-1],
        interfaces[-1],
        chemical,
        network,
        free,
        jump,
        times[-1],
        args.x_half,
        args.y_half,
        args.output_dir.name.replace("_", " ").title(),
    )
    make_diagnostics_figure(diagnostics_path, diagnostics)

    h = 2.0 * args.x_half / args.nx
    if not np.isclose(h, 2.0 * args.y_half / args.ny):
        raise ValueError("visualization requires equal x and y grid spacing")
    chemical_diag = diagnostics["chemical"]
    actin_diag = diagnostics["actin"]
    initial_actin_mass = float(diagnostics["initial_actin_mass"])
    actin_step_changes = np.diff(np.r_[initial_actin_mass, actin_diag[:, 2]])
    largest_actin_change_index = int(np.argmax(np.abs(actin_step_changes)))
    summary = {
        "grid": {
            "nx": args.nx,
            "ny": args.ny,
            "x_half": args.x_half,
            "y_half": args.y_half,
            "bulk_spacing": h,
        },
        "time": {
            "dt": args.dt,
            "steps": frames[-1] * args.output_every,
            "output_every": args.output_every,
            "final_time": float(times[-1]),
        },
        "interface": {
            "markers": args.nmarkers,
            "initial": interface_metrics(interfaces[0], h),
            "final": interface_metrics(interfaces[-1], h),
        },
        "fields": {
            "global_max_speed": float(np.hypot(u_all, v_all).max()),
            "final_pressure_min": float(p_all[-1].min()),
            "final_pressure_max": float(p_all[-1].max()),
            "final_chemical_min": float(chemical.min()),
            "final_chemical_max": float(chemical.max()),
            "final_network_actin_max": float(network.max()),
            "final_free_actin_max": float(free.max()),
            "final_jump_min": float(jump.min()),
            "final_jump_max": float(jump.max()),
        },
        "diagnostics": {
            "chemical_mass_first_reported": float(chemical_diag[0, 2]),
            "chemical_mass_final": float(chemical_diag[-1, 2]),
            "chemical_mass_change_from_first_report": float(
                chemical_diag[-1, 2] - chemical_diag[0, 2]
            ),
            "actin_mass_initial": initial_actin_mass,
            "actin_mass_final": float(actin_diag[-1, 2]),
            "actin_mass_relative_change": float(actin_diag[-1, 2] / initial_actin_mass - 1.0),
            "largest_absolute_one_step_actin_mass_change": float(actin_step_changes[largest_actin_change_index]),
            "largest_actin_mass_change_time": float(actin_diag[largest_actin_change_index, 1]),
            "gmres_converged_count": int(diagnostics["gmres_converged_count"]),
            "geometry_skip_warning_count": int(diagnostics["geometry_skip_warning_count"]),
            "mpi_bind_warning_count": int(diagnostics["mpi_bind_warning_count"]),
        },
        "artifacts": {
            "movie": movie_path.name,
            "final_figure": final_path.name,
            "diagnostics_figure": diagnostics_path.name,
        },
    }
    (args.output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
