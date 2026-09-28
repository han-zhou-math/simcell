# Cases for Figures 2, 3 and 4

One folder per figure. Each contains

- `inputs/<case>.par` — the exact `input.par` that the run used (copied from
  the run directory, not from a template);
- `cases.csv` — one row per case: its role in the figure, grid, marker
  count, box size, time step, completion state and last valid time, the
  build and SHA-256 of the executable that ran it, the location of the raw
  run in the study workspace, and the SHA-256 of the input;
- `figure/` — the accepted figure (PDF and PNG), its caption, and every
  table and audit behind the plotted numbers.

| Figure | Cases | Grid, markers, box, step | Builds used |
|---|---|---|---|
| 2 | 61 (7 + 9 + 6 per saturation level × 3 levels, plus field examples) | 128², 300, 80 µm, 5 s | `redesign_source128` (31), `execution_core` (30) |
| 3 | 14 (3 field/shape cases, strength response, width response) | 128², 300, 80 µm, 5 s | `redesign_source128` (11), `execution_core` (3) |
| 4 | 12 at 128² (neither / actin / pump / both across five pump profiles) + 4 at 256² (resolution check) | 128², 300, 100 µm, 5 s; 256², 300, 100 µm, 2.5 s | `figure4_source128` (12), `figure4_source256` (4) |

## Builds

All four builds compile the source in this repository (the tree was
synchronized with the frozen study copy on 27 September 2026). Only the
compile-time macros differ:

| Build | Macros | Executable SHA-256 (prefix) |
|---|---|---|
| `redesign_source128` | `-DSIMCELL_NRING=300` | `7cf46edaf0c33f1f` |
| `figure4_source128` | `-DSIMCELL_NRING=300 -DSIMCELL_X_HALF=1.25d0 -DSIMCELL_Y_HALF=1.25d0` | `2daf6231fc7139b7` |
| `figure4_source256` | `-DSIMCELL_NX=256 -DSIMCELL_NY=256 -DSIMCELL_NRING=300 -DSIMCELL_X_HALF=1.25d0 -DSIMCELL_Y_HALF=1.25d0` | `31fece10` |
| `execution_core` | `-DSIMCELL_NRING=300`, compiled from an earlier revision of the same source | `81ea2bcf2d2ac542` |

The `execution_core` runs predate the front/rear pump profile
(`dualchem_rear_width`, `dualchem_rear_amplitude_ratio`,
`dualchem_pump_start_time`) and the active actin-drive switch
(`stage14_active_actin_feedback`). Their inputs omit those namelist entries;
the current source gives them defaults that reproduce the earlier behaviour
exactly (a negative rear width selects the symmetric pump, and the drive is
on). The Figure 3 reuse audit confirmed that these were the only differences
between the reused inputs and new ones, and the current source runs those
inputs unchanged.

## Rerunning a case

```sh
make all SIMCELL_EXTRA_CPPFLAGS='<macros from the table>'
mkdir -p work/<case>/Data
cp cases/figureN/inputs/<case>.par work/<case>/input.par
cd work/<case> && OMP_NUM_THREADS=1 ../../imp > run.log 2>&1
```

Runs stopped before 2000 s where `cases.csv` says `stopped` or
`early_stop_boundary`; the figures use each case only up to its last valid
time, and the tables in `figure/` state the common cutoff used per panel.

## Raw output

Fields every 50 s and the per-step scalar logs for every case are in the
study workspace under the `source_run` path in `cases.csv`
(`/users/PBS0318/hzhou24/simcell/run_plan_discussion2/…/runs/<case>/`), about
1.3 GB in total. They are not versioned.
