# CUDA far-field grid resolution — Windows server reproduction (Option A, SPACING 1.0 → 0.5 Å)

- **Date / host**: 2026-09-05, Windows Server `117.50.175.42` (RDP 3389; SSH 22, user `administrator`). CPU host `10-60-245-8`, **NVIDIA RTX 5090**, driver 616.56. CUDA Toolkit v12.8 (nvcc 12.8.93) + v13.1 present; VS 2022 BuildTools (VC 14.44); rustc/cargo 1.98.1 (x86_64-pc-windows-msvc) in `C:\Users\Administrator\.cargo`.
- **Purpose**: apply the same far-field grid optimization (1 Å → 0.5 Å) to the **CUDA** engine and verify, by compiling + running on a real Windows/NVIDIA box, that the accuracy gain is inherited by the GPU far-field path.

## How the CUDA path is optimized
- In LKlight-GPU the 10–30 Å electrostatics far-field grid is **built once on the host by `src/grid_dna.rs` (`ReceptorField::build`)**; the CUDA kernels (`far_field.cu::far_field_kernel`, `full_score.cu::full_score_kernel`/`batch_full_score_kernel`) receive the grid (`phi`, `nx,ny,nz`, `origin`, `spacing`) and are **fully spacing-parametric** (trilinear gather). No kernel change is needed.
- The GPU repo's `src/grid_dna.rs` is byte-identical to the CPU grid repo (md5 `9b1c6397…`). So the same Option A patch applies:
  - `SPACING: f64 = 1.0 → 0.5`
  - generalise the six hardcoded `±32` scatter-window literals to `spread = ceil(FIELD_RMAX/spacing)+2` (fixes the latent truncation at non-1 Å spacing).
- Two CUDA builds produced on the Windows box (`cmd /c vcvars64.bat && cargo build --release --features cuda`):
  - `C:\lkwork\lkgpu-a\…\LKlight.exe` — Option A (0.5 Å), md5 `d5cb4bbd…`
  - `C:\lkwork\lkgpu-b\…\LKlight.exe` — baseline (1.0 Å), md5 `c4c1c720…`
  Both contain compiled kernels (`out/full_score.obj`, `out/farfield.lib`) and run `[gpu_score] CUDA full-pose scoring ACTIVE`.

## Validation (RTX 5090, CUDA ACTIVE)
Same 31-pose decoy set (seed 20260904) as the Mac/Linux CPU experiments; reference = deterministic exact (all-pairs) energies. Grid-limited scorers dna/pydock/cpydock.

| engine | dna max\|Δ\| | abs RMSD | mean\|Δ\| |
|---|---|---|---|
| GPU baseline (1.0 Å) | 9.601 | 4.094 | 3.379 |
| **GPU Option A (0.5 Å)** | **6.937** | **2.735** | **2.246** |
| (CPU grid 1.0 Å ref) | 9.600 | 4.094 | 3.379 |
| (CPU grid 0.5 Å ref) | 6.937 | 2.735 | 2.246 |

- The CUDA far field inherits the identical ~1.4× error cut (max|Δ| 9.6 → 6.94, RMSD 4.09 → 2.74).
- **pydock / cpydock**: Windows-GPU energies are **bit-identical** to the Linux-CPU grid at the same spacing (max diff 0.000).
- **dna**: Windows-GPU vs Linux-CPU differ by at most **6.1e-4** (mean 7e-5) — the documented f32 near-field rounding of the CUDA batch path; far-field behaviour matches.

## Raw outputs
- `win_gpu_A_perpose.csv`, `win_gpu_B_perpose.csv` (this directory; also `C:\lkwork\` on the server).
- Sources + binaries retained in `C:\lkwork\lkgpu-a`, `C:\lkwork\lkgpu-b` on the server; build helper `C:\lkwork\lk_build_cuda.bat`; pose list `poses.tsv`; driver `run_win_poses.ps1`.

## Note
This validates that the CUDA engine *can* carry the 0.5 Å optimization with the same accuracy as the CPU grid. Whether to adopt 0.5 Å as the shipped default (it raises the per-setup field-build cost ~8×; steady-state per-pose docking is unaffected because the field is built once) is a release decision pending user confirmation.
