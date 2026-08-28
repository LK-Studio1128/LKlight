# LKlight v1.1.0

**DFIRE2 non-protein-ligand fix + multi-scenario test hardening.** This release
completes the multi-module / multi-scenario validation campaign and fixes the
last known crash class in the scoring layer.

## Fixed

- **DFIRE2 panic on non-protein ligands** (`src/dfire2.rs`): `score rec lig dfire2`
  panicked with `index out of bounds` when the ligand contained no atoms in the
  DFIRE2 atom-type dictionary (e.g. DNA ligands in 1AZP / p53-DNA 1DIZ).
  The ligand residue-offset computation indexed `[0]` on an empty vector; it now
  uses safe `match (last, first)` destructuring with graceful degradation.
  Regression-verified: 1AZP dfire2 = −83.70; full p53-DNA GSO run completes.
- Carried hardening from the BM5 equivalence work: robust atom typing and ANM
  stride guards (`cpydock`, `ddna`, `mj3h`, `sd`, `sipper`), PDB free-text
  metadata tolerance (Fixes 5–6).

## Tested (new multi-scenario suite)

- **12 scoring functions × 6 biological complexes** (protein–DNA 1AZP / 1DIZ
  p53–DNA, protein–protein 2OOB, Ab–lysozyme 1VFB, Ab–peptide 1DQJ,
  viral–host 6M0J RBD–hACE2): 78+ score combinations, zero crashes, zero
  timeouts.
- **Parameter sweep**: wall-clock scales linearly with glowworms (25→400) and
  GSO steps (10→200); convergence by 50–70 steps — 100 steps remains the
  practical default.
- **All 16 CLI subcommands** exercised end-to-end.
- 29/29 unit tests pass.

## Downloads

| Platform | File | Binary |
|---|---|---|
| macOS Apple Silicon | `LKlight-v1.1.0-mac-arm64.zip` | aarch64, native |
| Linux x86-64 | `LKlight-v1.1.0-linux-x86_64.tar.gz` | musl **static-PIE**, runs on any distro |
| Windows x64 | `LKlight-v1.1.0-win-x64.zip` | PE32+ console |

Each archive contains the binary plus `README.md`, `LICENSE` (GPL-3.0),
`NOTICE` and `CHANGELOG.md`.

## Citation & DOI

This release is archived on Zenodo; the DOI badge in the README will update
automatically once Zenodo registers the new version DOI.

LKlight is a derivative work of the LightDock `lightdock-rust` Rust baseline
(GPL-3.0). See `NOTICE` for full attribution.
