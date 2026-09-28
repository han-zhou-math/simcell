# SimCELL dev-minimal

SimCELL simulates the migration of a single cell in two dimensions. The cell
is a closed, elastic, water-permeable membrane immersed in Stokes fluid.
Inside it, F-actin and G-actin are transported by the cytosol and exchanged
by polymerization at the membrane and depolymerization in the bulk. A solute
diffuses on both sides of the membrane, leaks through it, and is actively
pumped across it. The resulting osmotic water flux and actin stresses move
the membrane, and the cell migrates along a periodic channel. The numerical
method is an implicit immersed boundary method for the fluid–membrane
problem coupled to a Cartesian-grid advection–diffusion solver for moving
domains; a summary is given at the end of this file.

The code is a compact Fortran/C++ program: 33 Fortran solver files (28 in
the root, 5 in `dualchem/`), a 4-file C++ multigrid backend, and 10 unit
tests, built into one serial executable `imp`. This is the exact source that
produced Figures 2–4 of the osmosis–actin manuscript; the cases, inputs and
accepted figures are in `cases/`. Readers who want to run a
case should go to **Build** and **Run**; readers who want the numbers behind
a figure should go to **Inspecting results**.

## Contents

| Path | What it is |
|---|---|
| `fmain.f90` | Program entry point; one time step is fluid–membrane, then solute, then actin |
| `fsisolve.f90`, `linsys.f90`, `IBmod.f90`, `IBforce.f90`, `brinkman_solver_mod.f90` | Implicit immersed boundary fluid–membrane solver (FFT + banded Stokes solve, PETSc GMRES) |
| `dualchem/` | Solute advection–diffusion on the moving interior/exterior domains; `AdvecDiff2d-semiperiodic/` is the C++ multigrid backend |
| `actin_*.f90`, `fsi_actin_*_mod.f90` | F/G-actin transport, polymerization law, and actin feedback on the fluid and membrane |
| `chemical_*_mod.f90`, `osmotic_feedback_mod.f90`, `fsi_adhesion_mod.f90`, `fsi_external_load_mod.f90` | Membrane flux laws, osmotic slip, adhesion, and external loads |
| `parameters.f90` | Compile-time grid, marker count and box size; runtime namelist |
| `input.par` | Default 100-step demonstration input (namelist `PARAM`) |
| `cases/figure2/`, `cases/figure3/`, `cases/figure4/` | For each figure: the exact inputs of every case, a case table with build and run records, and the accepted figure with caption and tables |
| `tools/visualize_results.py` | Plots the frames and log of one run |

Only these files are versioned. Unit tests, convergence drivers, notes and
papers exist in the working copy but are not part of the repository.

Generated output goes to `Data/` (current run) and `output/` (archived
studies); both are ignored by git.

## Requirements

- PETSc (3.19 or later; 3.19.6 is what the study used) built with the same
  Fortran compiler you will use here. The Makefile reads compilers and flags
  from PETSc's `pkg-config` file.
- FFTW 3, BLAS and LAPACK.
- gfortran, g++ (C++11), GNU make.
- Python 3 with NumPy and Matplotlib for the visualization and figure
  scripts; `ffmpeg` if you want movies.

On the OSC cluster the study used a prepared environment:

```sh
source /users/PBS0318/hzhou24/simcell/.deps/activate-simcell.sh
```

which loads `gcc/12.3.0` and `openmpi/5.0.2` and exports `PETSC_DIR` and
`PETSC_ARCH`. Elsewhere, set those two variables yourself:

```sh
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-opt
```

## Build

```sh
make print   # show the compilers and flags resolved through PETSc
make all     # clean, compile and link ./imp
make imp     # incremental rebuild
make clean   # remove objects, modules, the executable and test binaries
```

The grid, marker count and box size are fixed at compile time. Defaults are
128 × 128 cells, 200 membrane markers, and a box of half-width 1 in code
units (80 µm with the 40 µm length scale). All figure runs used 300 markers,
so override the defaults with preprocessor macros:

```sh
# Figures 2 and 3: 128² cells, 300 markers, 80 µm box
make all SIMCELL_EXTRA_CPPFLAGS='-DSIMCELL_NRING=300'

# Figure 4: 128² cells, 300 markers, 100 µm box
make all SIMCELL_EXTRA_CPPFLAGS='-DSIMCELL_NRING=300 -DSIMCELL_X_HALF=1.25d0 -DSIMCELL_Y_HALF=1.25d0'

# Figure 4 resolution check: 256² cells, 300 markers, 100 µm box
make all SIMCELL_EXTRA_CPPFLAGS='-DSIMCELL_NX=256 -DSIMCELL_NY=256 -DSIMCELL_NRING=300 -DSIMCELL_X_HALF=1.25d0 -DSIMCELL_Y_HALF=1.25d0'
```

| Macro | Meaning | Default |
|---|---|---|
| `SIMCELL_NX`, `SIMCELL_NY` | Cells per direction (powers of two; the multigrid requires it) | 128 |
| `SIMCELL_NRING` | Number of membrane markers | 200 |
| `SIMCELL_X_HALF`, `SIMCELL_Y_HALF` | Half-width of the box in code length units | `1.0d0` |

A change of macros requires `make all` (a full rebuild), not `make imp`.

## Run

`imp` reads `input.par` from the current directory and writes into `./Data/`,
which must exist. Run each case in its own directory:

```sh
mkdir -p work/demo/Data
cp input.par work/demo/
cd work/demo && ../../imp | tee run.log
```

Set `OMP_NUM_THREADS=1` (and the BLAS equivalents) — the solver is serial and
oversubscribed threads slow it down. A 128² step takes about 5 s on one core;
a 256² step about 14 s.

### Input file

`input.par` is a Fortran namelist. All coefficients are **dimensionless**;
the `scale_*` entries only record the physical scales used for
postprocessing. The conversion from physical parameters is documented in
the nondimensionalization note kept with the manuscript. The most important
entries are:

| Entry | Meaning |
|---|---|
| `dlt`, `ntmax`, `nfreq` | Time step (code units), number of steps, and field-output interval in steps. `dlt = 0.005` with `scale_time_s = 1000` is 5 s |
| `nu` | Dimensionless viscosity |
| `clstiff`, `membrane_reference_metric` | Membrane tension coefficient and reference stretch (0 = the zero-reference-length law) |
| `initial_cell_radius`, `initial_cell_axis_ratio` | Initial shape in code length units |
| `stage12_actin_eta`, `stage12_actin_eta_s`, `stage12_actin_k_sigma`, `stage12_actin_dc`, `stage12_actin_gamma` | Cytosol–network drag, focal adhesion, F-actin stress, G-actin diffusivity, depolymerization rate |
| `stage12_actin_jc`, `stage12_actin_theta0`, `stage12_actin_dw`, `stage12_localized_polymerization` | Polymerization flux amplitude, saturation concentration, leading-edge width, and the switch for the localized saturating law |
| `dualchem_diffusion`, `dualchem_kc`, `dualchem_kp`, `dualchem_initial_concentration` | Solute diffusivity, passive leak, pump strength, and initial concentration |
| `enable_osmotic_feedback`, `water_osmotic_mobility`, `water_stress_mobility` | Osmotic slip switch and the two hydraulic mobilities (osmotic jump, membrane traction) |
| `stage14_adhesion`, `stage14_external_load_x` | Membrane–substrate resistance and an optional external load |
| `stage14_recenter_interface` | Comoving frame: remove the mean x-translation of the membrane each step (the x direction is periodic) |
| `scale_length_um`, `scale_time_s`, `scale_velocity_um_s`, `scale_stress_pa`, `scale_concentration_millimolar` | Physical scales recorded for postprocessing |

The inputs used for the paper are in `cases/figures_2_3_4/*/inputs/*.par`
and are the best starting point for a new case.

### Running a published case by hand

```sh
# 1. build with the macros the campaign used (see its README / build.sh)
make all SIMCELL_EXTRA_CPPFLAGS='-DSIMCELL_NX=128 -DSIMCELL_NY=128 -DSIMCELL_NRING=300 -DSIMCELL_X_HALF=1.25d0 -DSIMCELL_Y_HALF=1.25d0'

# 2. copy the case input and run
mkdir -p work/F4_L3_BASE/Data
cp cases/figure4/inputs/F4_L3_BASE.par work/F4_L3_BASE/input.par
cd work/F4_L3_BASE && OMP_NUM_THREADS=1 ../../imp > run.log 2>&1
```

`cases/figureN/cases.csv` lists every case of that figure with its role in
the figure, grid, marker count, box, time step, completion state, last valid
time, the build that ran it, and the executable and input hashes; the
matching compile macros are in `cases/README.md`. The inputs are the
`input.par` files copied from the actual run directories. Runs that
`cases.csv` marks `stopped` ended before 2000 s on a solver failure or a
boundary encounter; the figures use those cases only up to their last valid
time.

## Output

**Standard output** carries one diagnostic line per step for each stage.
Lines starting with these tags are the ones to watch:

| Tag | Content |
|---|---|
| `CAMPAIGN_RUNTIME`, `CAMPAIGN_COMMITTED` | Grid, fixed-step flag; accepted step, time, and minimum concentrations |
| `STAGE07_DUALCHEM_MASS`, `STAGE07_PHYSICAL_JUMP` | Total solute and the membrane concentration jump |
| `STAGE12_ACTIN_MASS` | Total F + G actin in the cell |
| `STAGE13_BRINKMAN` | Fluid solve counts and inner iterations |
| `STAGE14_ACTIN_FEEDBACK` | Maximum actin drag, body force and membrane stress |
| `STAGE15_COMOVING_SHIFT` | Per-step and cumulative frame shift when recentering is on (the cumulative value is the lab-frame displacement) |
| `F4_STARTUP` | Every step, in physical units: time (s); lab-frame centroid x and y (µm, the x value includes the comoving shift); cell area; total F-actin, total G-actin and the F fraction; front and rear G-actin and their polarity; front and rear F-actin and their polarity; pump start time |
| `REDESIGN_ACTIVE_FORCE` | Whether the active actin drive is on and the maximum active force terms (used by the Figure 4 controls) |
| `ERROR STOP …`, `CAMPAIGN_FIXED_STEP_REJECT` | A solver failure, or the membrane moving more than one cell in a step; the run stops |

**`Data/`** receives a frame every `nfreq` steps, numbered `0000, 0001, …`:

| File | Content | Format |
|---|---|---|
| `frun.ib.NNNN` | Marker coordinates | text, `x y` per marker |
| `frun.u.NNNN`, `frun.v.NNNN`, `frun.p.NNNN` | Velocity components and pressure on the MAC grid | raw `float64`, Fortran order |
| `frun.c.NNNN` | Solute concentration (both sides, cell centers) | raw `float64` |
| `frun.n.NNNN`, `frun.g.NNNN` | F-actin and G-actin (cell centers, meaningful inside the cell) | raw `float64` |
| `stage07.chemical.final.bin`, `stage07.jump.final.bin`, `stage12.network.final.bin`, `stage12.free.final.bin` | Final solute field and jump, final F- and G-actin | raw `float64` |

Read a field in Python with
`np.fromfile(path, dtype=np.float64).reshape((nx, ny), order="F")`
(`nx × (ny−1)` for `v`). Convert to physical units with the `scale_*`
entries of the input.

## Inspecting results

**Quick look at one run.** `tools/visualize_results.py` reads a run's
`Data/` frames and log and writes field maps, membrane overlays, and time
series (and a movie if `ffmpeg` is available):

```sh
python3 tools/visualize_results.py --data-dir work/demo/Data --log work/demo/run.log \
    --output-dir work/demo/figures --nx 128 --ny 128 --nmarkers 200 --dt 0.005 \
    --output-every 5 --x-half 1.0 --y-half 1.0
```

Pass the values that match the build and input (`--nmarkers 300` for any
figure case; add `--x-half 1.25 --y-half 1.25` for Figure 4; `--nx 256 --ny
256 --dt 0.0025` for the resolution check).

**Accepted figures and tables.** `cases/figureN/figure/` holds the accepted
figure with its caption and the tables behind every plotted number:

| Figure | Figure files | Tables / audits |
|---|---|---|
| 2 | `Figure_2.pdf`, `Figure_2.png`, `Figure_2_caption.md` | selected cases, velocity source, exclusions, conservation, `conservation_audit.json` |
| 3 | `Figure_3.pdf`, `Figure_3_support.pdf`, `Figure_3_caption.md` | cases, strength response, width response, displacement, conservation, permeability audit, per-case mass histories |
| 4 | `Figure_4.pdf`, `caption.md` (plus the interaction-preview variant) | interaction summary and time series, 256² mesh check |

`cases/README.md` explains the case tables, the four builds, and where the
raw run output is stored.

**Conventions used in the figures.** Speeds are in µm/h over a stated time
interval, pump strengths in µm/s, concentrations in mM. Cases are compared on
a common valid time interval, never on each case's last available frame.
Total solute and total actin are reconstructed from the saved fields with
cut-cell area fractions and one-sided values; the F-actin fraction is not a
conservation measure.

## Checking a build

After building, run the default `input.par` (100 steps, about ten minutes on
one core) and confirm that `STAGE07_DUALCHEM_MASS` and `STAGE12_ACTIN_MASS`
stay constant to solver tolerance, that every `CAMPAIGN_COMMITTED` line
reports nonnegative minimum concentrations, and that no `ERROR STOP` occurs.
The Makefile also lists `check_*` targets for unit tests whose sources are
not part of this repository.

## Numerical method

Fluid and membrane: Stokes–Brinkman on a MAC grid, membrane markers with the
four-point immersed boundary kernel, implicit midpoint update of the membrane
solved by an approximate-Newton/GMRES iteration (Mori & Peskin 2008;
Yao & Li, J. Theor. Biol. 2025). Solute and actin: the moving-domain problem
is extended to the fixed box with an unknown jump in normal derivative on the
membrane; backward Euler with centered conservative differences, a quadratic
correction function at cut cells, geometric multigrid for the bulk and
matrix-free GMRES for the interface density (Zhou, Mori & Yao 2026). The
full description with equations is the manuscript's numerical method
section.

## Provenance of the figure results

The Figure 2–4 runs used four executables, all compiled from this source
(the versioned tree was synchronized on 27 September 2026 with the frozen
study copy `run_plan_discussion2/redesign_2026_09_18/source128`):

| Build | Used by | Compile macros | Executable SHA-256 (prefix) |
|---|---|---|---|
| `redesign_source128` | 31 Figure 2 cases, 11 Figure 3 cases | `NRING=300` | `7cf46edaf0c33f1f` |
| `execution_core` | 30 Figure 2 cases, 3 Figure 3 cases | `NRING=300`; an earlier revision of the same source, without the front/rear pump and active-drive inputs | `81ea2bcf2d2ac542` |
| `figure4_source128` | 12 Figure 4 cases | `NRING=300`, `X_HALF=Y_HALF=1.25d0` | `2daf6231fc7139b7` |
| `figure4_source256` | 4 Figure 4 resolution-check cases | `NX=NY=256`, `NRING=300`, `X_HALF=Y_HALF=1.25d0` | `31fece10` |

`cases/figureN/cases.csv` gives the build and full hash for every case. The
current source runs the `execution_core` inputs unchanged: the namelist
entries they lack take defaults that reproduce the earlier behaviour (a
negative `dualchem_rear_width` selects the symmetric pump exactly, and the
active actin drive is on). A fresh build will not reproduce the hashes bit
for bit (build paths and timestamps differ) but is the same source.

## History

Imported from the Stage 15 worktree (`adbc29c`) as a minimal functional
core, followed by: deterministic ghost handling; the moving-interface Robin
condition with relative normal velocity and true-residual GMRES acceptance;
one explicit common solute diffusivity; the hybrid chemical GMRES; osmosis
studies, external loads and the fixed interface stencil; the coupled
osmosis–actin study; and the front/rear pump profile, active-drive switch,
initial actin concentrations, fixed-step enforcement and campaign
diagnostics added for the Figure 2–4 campaigns. Only the Figure 2–4 cases
are versioned under `cases/`; earlier study folders remain on disk
untracked.

Sources carry SPDX BSD-3-Clause headers.
