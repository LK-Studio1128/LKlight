# Far-field grid resolution — Linux server reproduction (Option A, SPACING 1.0 → 0.5 Å)

- **Date / host**: 2026-09-05, Linux server `117.50.160.76` (Ubuntu 22.04 / RTX 3080 Ti 12 GB / CUDA 13.2 / rustc 1.98.0, 12 cores). Logged in as `ubuntu`.
- **Purpose**: independently re-verify on a real Linux x86 server the macOS (Apple Silicon) result that halving the far-field grid spacing (1 Å → 0.5 Å) cuts the grid-vs-exact absolute error ~1.4×, to confirm no platform-specific behaviour before any paper statement.
- **Sources**:
  - Exact (all-pairs reference): `~/LKlight` repo, `cargo build --release` (27.7 s), bin `~/LKlight/target/release/LKlight`.
  - Grid baseline: `~/lk_acc/lklight-grid-src` = GitHub `LK-Studio1128/LKlight-grid` @ HEAD `1f3bfb0` (`grid_dna.rs` md5 `9b1c6397f787392eabad1d114f548aa7`, identical to shipped + to Mac working tree). Built `cargo build --release` (CPU, no cuda feature) → `~/lk_acc/LKlight-grid-linux-BASE`.
  - Grid Option A: copy of the above with `grid_dna.rs` patched → `~/lk_acc/LKlight-grid-linux-A`. Patch = `SPACING 1.0→0.5` + generalise the six hardcoded `±32` scatter-window literals to `spread = ceil(FIELD_RMAX/s)+2` (the latent bug that truncates the far field at non-1 Å spacing). Backup `grid_dna.rs.ORIG`/`grid_dna.rs.PINNED` retained.
- **Harness**: `equiv_harness_server.py` (env-driven clone of the Mac `equiv_harness.py`, same deterministic decoy seed 20260904), 1 native + 30 decoys = 31 poses, 12 scoring functions. Grid binary compared against the exact reference.
- **Raw outputs**: `equiv_linux_BASE_{summary,perpose}.csv`, `equiv_linux_A_{summary,perpose}.csv` (this directory), also mirrored to `/home/ubuntu/lk_acc/` on the server.

## Results (absolute far-field error, energy units, 31 poses)

| metric | Linux BASE (1 Å) | Linux A (0.5 Å) | Mac A (0.5 Å) |
|---|---|---|---|
| max \|Δ\| | 9.600 | **6.937** | 6.94 |
| abs RMSD | 4.094 | **2.735** | 2.74 |
| mean \|Δ\| | 3.379 | **2.246** | — |
| 9 non-far-field fn | bit-identical | bit-identical | bit-identical |
| dna / pydock Spearman | 1.000 | 1.000 | 1.000 |
| cpydock Spearman | 0.9770 | 0.9879 | 0.988 |

**Cross-platform determinism**: for every pose, the exact and grid energies (both BASE and A) are **bit-identical between macOS arm64 and Linux x86** (max platform difference 0.00e+00 across dna/pydock/cpydock, and the A-vs-B improvement per pose is likewise identical). The CPU scoring + far-field grid path is fully deterministic across the two platforms.

**Conclusion**: Option A's accuracy gain is reproduced exactly on the real Linux server (max|Δ| 9.60 → 6.94, absRMSD 4.09 → 2.74). The 1 Å grid "reference-level accuracy" claim and the ~1.4× improvement figure are platform-independent and now carry a genuine server-measured datapoint.

*Note: Option A is a code-level experiment only and is NOT being adopted into the engine or paper; the shipped 1 Å grid behaviour is unchanged. Source trees on server and Mac were restored/left unmodified after the run.*
