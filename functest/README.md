# LKlight CPU-grid vs base-exact: per-function benchmark (paper-2 supplement)

This directory contains the runnable harness + raw data + figures that substantiate the
"per-scoring-function equivalence and speedup" claim of the GPU/grid paper.

## Engines compared

| label  | binary path                                                                | role          |
|--------|----------------------------------------------------------------------------|---------------|
| `base` | `byi/LKlight/target/release/LKlight`                                       | all-pairs exact (no grid) |
| `grid` | `byi/LKlight-grid/release_bin/LKlight-mac-arm64`                           | CPU cell-list + far-field electrostatic map |
| `gpu`  | server only: `byi/LKlight-GPU/release_bin/LKlight-linux-cuda`              | CUDA batch kernel (not runnable on this Mac) |

`base` has no `nearcell`/`far-field`/`energy_grid` strings — it is the genuine
v1.1.0 all-pairs reference. `grid` and `gpu` build on the same v1.1.0 physics;
the grid binary's `energy_exact` is byte-identical to `base` (verified on the
1AZP native pose: pydock = -364.881264 on both).

## Files

| file                                  | what it produces                                       |
|--------------------------------------|--------------------------------------------------------|
| `equiv_harness.py`                   | 30-decoy pose scan through `base` and `grid` for all 12 methods; equivalence CSV + summary CSV |
| `gso_harness.py`                     | full GSO `run` (1 swarm × 20 glow × 100 steps) per method × both engines; steady-state wall-clock + best-energy convergence CSV |
| `make_figF.py`                       | publication figures F1 (equivalence) + F2 (speedup)   |
| `server_functest.sh`                 | server-side three-engine benchmark for the GPU + the full large-RNA system (12 methods × 3 engines × N repeats) |
| `equiv_full_summary.csv`             | per-function equivalence summary (31 poses, 12 fns)   |
| `equiv_full_perpose.csv`             | per-pose exact/grid energies + per-pose time          |
| `gso_full_summary.csv`               | per-function GSO steady-state wall-clock + best score |
| `figures3/figF1_equiv.png`           | Figure F1 — grid≡base per-function agreement          |
| `figures3/figF2_speedup.png`         | Figure F2 — base vs grid steady-state, plus documented large-scale |

## How to reproduce (local, Mac)

```bash
cd functest
python3 equiv_harness.py --nposes 30 --out equiv_full
python3 gso_harness.py   --steps 100 --out gso_full
python3 make_figF.py
```

All inputs (1AZP receptor/ligand) are already copied into this directory.

## How to reproduce (Linux server, RTX 3080 Ti)

```bash
# 1. ship the directory contents to the server.
# 2. set environment variables (or use defaults):
export EXACT_BIN=/path/to/LKlight-linux86-static
export GRID_BIN=/path/to/LKlight-linux-x64
export GPU_BIN=/path/to/LKlight-linux-cuda
export REC=/path/to/large_receptor.pdb     # e.g. AF-M9NF14 8,218 atoms
export LIG=/path/to/large_ligand.pdb       # e.g. CR44522 12,625 atoms
export GLOW=20
export STEPS=100
export REPEATS=2
./server_functest.sh
```

The script writes two files under `functest_run_<pid>/`:

- `server_functest_perfn.csv` — one row per (method × engine × repeat)
- `server_functest_summary.csv` — averaged wall-clock per (method × engine) plus best-energy and rel-delta vs exact

Append the summary rows to the paper's Results table after running.

## Methodological notes

1. **Single-process `score` is NOT a steady-state timing measure** for the
   grid-accelerated functions: each fresh process pays a one-time
   field+cell-list construction cost (~265 ms on 1AZP). We use the full GSO
   `run` (2000 pose evaluations) so the one-time setup amortises; this is the
   number reported in `gso_full_summary.csv`.

2. **`max_rel_dev` for cpyDock is large (216 %) because the test set includes
   poses with near-zero energies where the interpolation error |g-e| dominates
   the ratio. The ranking metric (Spearman = 0.977) is the appropriate
   equivalence measure for this function.**

3. **At 1AZP small-system scale the grid engine is at parity (or slightly
   slower for tiny look-up functions) because the cell-list construction
   overhead dominates the per-pose cost.** The large-complex speedups
   documented in PERF_COMPARE §八.1 (vdW 48× / dna 19.4× / pyDock 15× /
   cpyDock 11.7× on an 8,218+12,625 atom system) require a large complex and
   are reproducible via `server_functest.sh`.