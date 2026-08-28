# Changelog

All notable changes to LKlight are documented in this file.
LKlight is a derivative work of the LightDock `lightdock-rust` Rust baseline.
Changes listed below are relative to that upstream baseline.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [1.1.0] — 2026-08-29

### Fixed (found during multi-scenario testing, 2026-08-29)

- **[Fix 7] DFIRE2 panic on non-protein ligands** (`src/dfire2.rs`):
  `score rec lig dfire2` panicked with `index out of bounds: the len is 0`
  whenever the ligand contained no atoms in the DFIRE2 atom-type dictionary
  (e.g. DNA ligands in 1AZP/1DIZ — DFIRE2 ships only the 20 standard amino
  acids): the ligand residue-offset computation indexed `[0]` on an empty
  vector. Now uses safe `match (last, first)` destructuring and degrades
  gracefully (offset 0). Regression-verified: 1AZP dfire2 = -83.70, p53-DNA
  full GSO run completes normally.

- Carried (previously uncommitted) hardening from the BM5 equivalence work:
  robust atom typing and ANM stride guards in `src/cpydock.rs`, `src/ddna.rs`,
  `src/mj3h.rs`, `src/sd.rs`, `src/sipper.rs`, and PDB free-text metadata
  tolerance in `src/bin/lightdock.rs` (Fixes 5–6 below).

### Tested (multi-module / multi-scenario suite, 2026-08-29)

- 12 scoring functions × 6 biological complexes (protein-DNA 1AZP/1DIZ p53-DNA,
  protein-protein 2OOB, Ab-antigen 1VFB / Ab-peptide 1DQJ, viral-host 6M0J
  RBD-hACE2): 78+ score combinations, no timeout, no panic.
- Parameter sweep: runtime scales linearly with glowworms (25→400) and steps
  (10→200); GSO converges by 50–70 steps; 100 steps is the practical default.
- All 16 CLI subcommands exercised end-to-end (setup/run/rank/cluster/top/
  generate/trajectory/map_contacts/diameter/gso_to_csv/reference_points/
  pipeline/...).

### Fixed (found during BM5 equivalence benchmarking)

- **[Fix 5] PDB parsing of free-text metadata records** (`src/bin/lightdock.rs`):
  `pdbtbx 0.11` raises an invalidating error on free-text metadata lines such as
  `REMARK DATE:23-Dec-2018` (it attempts to parse them as `usize`), which made
  LKlight crash on most real-world PDB files downloaded from the RCSB or produced
  by `reduce`. `open_pdb_padded` now strips metadata records (REMARK/USER/HEADER/
  TITLE/COMPND/SOURCE/KEYWDS/EXPDTA/AUTHOR/REVDAT/JRNL/FORMUL/HET) and retries the
  parse before failing; atom coordinates are never modified.

- **[Fix 6] Swarm starting positions too far from the receptor** (`src/bin/lightdock.rs`):
  Initial ligand poses were placed on a sphere of `max_atom_radius + 10 Å`, which
  for typical complexes put the ligand ~40–45 Å from the receptor centre — too far
  for blind docking to converge (success rate 0 in BM5 tests). The sphere radius is
  now `mean_atom_radius + swarm_radius` (default `swarm_radius` reduced 10 → 3 Å),
  placing initial poses near the molecular surface, matching the official
  LightDock reference-point convention (~15–25 Å). Use `--swarm-radius` to adjust.

---

## [1.0.0] — 2025

### Critical Bug Fixes (relative to lightdock-rust baseline)

- **[Fix 1] DFIRE/DFIRE2/dDNA parameter embedding** (`src/dfire.rs`, `src/dfire2.rs`, `src/ddna.rs`):
  The upstream baseline loaded DFIRE parameter matrices from external files at
  runtime. When those files were absent (as in any standalone binary distribution),
  the program panicked with `Unable to open DFIRE parameters: NotFound`. LKlight
  embeds all parameter matrices at compile time using Rust constant arrays,
  eliminating any external file dependency and making DFIRE/DFIRE2/dDNA scoring
  fully portable for the first time.

- **[Fix 2] ANM stride computation** (`src/pisa.rs`, `src/ddna.rs`, `src/cpydock.rs`, `src/sd.rs`):
  The upstream code used a fixed ANM atom stride that did not match the actual
  1D storage layout of mode vectors (`nmodes.len() / (3 * n_modes)`), causing
  out-of-bounds reads and incorrect coordinate displacements in all ANM-enabled
  docking runs. Corrected to the proper stride and added a bounds guard
  `if i_atom >= nm_n { break; }` throughout.

- **[Fix 3] Non-standard residue handling** (`src/dfire.rs`):
  DFIRE scoring panicked on any residue outside the standard 20 amino acids
  (e.g. ligands, modified residues, DNA nucleotides). LKlight returns a penalty
  value of 999 for unknown residue types and emits a warning instead of crashing.

- **[Fix 4] ANM atom-count assertion** (`src/simulator.rs`):
  Removed `assert_eq!(model.atom_count(), anm_atoms)`. ANM in LightDock covers
  only Cα atoms, not all heavy atoms; the assertion was always incorrect and
  caused every ANM-enabled run to abort immediately.

### Added: Performance Optimizations

- **[H1] rayon outer-loop parallelization** (`src/pydock.rs`, `src/dna.rs`,
  `src/cpydock.rs`, `src/dfire.rs`, `src/dfire2.rs`, `src/sd.rs`):
  Parallelized the receptor-atom outer loop using `rayon::par_iter()`. Energy
  accumulation (ELEC + VDW / DFIRE bins) is computed in parallel; the
  interface-flag phase (INTERFACE_CUTOFF = 3.9 Å, <5% of runtime) is kept
  sequential to avoid synchronization overhead. Benchmark improvement:
  pydock 7693 ms → 290 ms (26.5×), dfire 935 ms → 33 ms (35×).

- **[H2] SIMD-friendly hot paths**:
  Refactored hot loops and data layout so the LLVM backend can auto-vectorize
  distance calculation and energy accumulation when appropriate. Published
  binaries use portable CPU baselines; local benchmark builds may additionally
  enable `target-cpu=native`.

- **[F1] sqrt_vdw_charges precomputation** (`src/pydock.rs`, `src/cpydock.rs`,
  `src/sd.rs`, `src/dna.rs`): Precomputed `vdw_charge.sqrt()` at initialization
  time; hot-path `(rec.vdw_charges[i] * lig.vdw_charges[j]).sqrt()` replaced
  with `rec.sqrt_vdw_charges[i] * lig.sqrt_vdw_charges[j]`.

- **[F2] SD spatial hash grid** (`src/sd.rs`): Implemented a `thread_local!`
  `HashMap<(i32,i32,i32), Vec<usize>>` 3D grid (cell = 9 Å) for the SD scoring
  function. With a 9 Å cutoff much shorter than typical protein diameters, the
  grid provides genuine O(N²)→O(N) sparsification querying only 27 neighboring
  cells per receptor atom.

- **[F3] BufWriter I/O batching** (`src/swarm.rs`): Replaced per-line
  `write!()` syscalls in `Swarm::save()` with `BufWriter<File>`, reducing
  filesystem syscall count by ~100× per GSO output file.

- **[F4] Stack-allocated quaternion rotate** (`src/qt.rs`): Changed
  `rotate()` return type from `Vec<f64>` (heap, 1 allocation per call)
  to `[f64; 3]` (stack), eliminating allocations in the GSO rotation hot path.

- **[F5] Glowworm stack translation** (`src/glowworm.rs`): Changed
  `Glowworm::translation` from `Vec<f64>` to `[f64; 3]`.

- **[G1] Swarm scratch buffers** (`src/swarm.rs`): Added `pos_scratch:
  Vec<[f64; 3]>` and `rot_scratch: Vec<Quaternion>` fields to `Swarm`,
  reused per `movement_phase()` call via `resize + fill`.

- **[G2] GSO movement parallelization** (`src/swarm.rs`): Parallelized
  `movement_phase()` with `rayon::par_iter_mut()`. Random numbers are
  pre-generated before the parallel region to preserve determinism.

- **thread_local! scoring scratch buffers**: All scoring functions use
  `thread_local! { static SCRATCH: RefCell<(Vec<f64>, Vec<u8>)> }` for
  coordinate and interface flag buffers, eliminating per-call heap allocations.

### Changed: Spatial Grid Retraction (G3/G4 → H1)

- Initial attempts to apply 10 Å/±3-cell `HashMap` spatial grids to
  `pydock.rs`, `dna.rs`, `cpydock.rs` (G3/G4) produced **performance
  regression**: with 30 Å electrostatic cutoffs covering most of a small
  protein ligand, 343 HashMap queries per receptor atom (7³ cells) far
  exceeded the cost of a direct O(N²) scan (~221 ligand atoms × 2 ns vs.
  343 × 50 ns). These grids were retracted; the final H1 approach uses
  `rayon` parallelization instead.

- Similarly, DFIRE's initial 15 Å grid (27-cell, ±1) was retracted (I1/I2)
  for the same reason: the 15 Å cutoff is comparable to a small protein,
  so the grid prunes almost nothing while imposing HashMap overhead.
  The sd.rs 9 Å grid (F2) is retained — its short cutoff yields genuine
  sparsification.

### Added: Branding and Interface

- **Renamed** project from `lightdock-rust` to `LKlight`; binary renamed
  from `lightdock` to `LKlight`.
- **Unified CLI** entry point in `src/bin/lightdock.rs` with subcommands:
  `setup`, `run`, `generate`, `cluster`, `rank`, `top`, `filter`, `score`,
  `trajectory`, `pipeline`.
- Added `NOTICE`, `CHANGELOG.md`, `CONTRIBUTING.md` for GPL compliance.

### Testing

- 29 unit tests (`cargo test --lib`): **29/29 pass**
- 160 numerical integration tests vs. Python LightDock reference: **160/160 pass**
- Floating-point equality assertions updated to approximate tolerance
  (ε = 10⁻⁸) to accommodate accumulation-order differences from parallel
  reduction.
