# LKlight v1.2.1

**Nucleic-acid scoring safety patch.** This patch prevents protein–DNA / RNA
jobs from silently completing with an unsupported scoring function and a
constant zero / `-0` score.

## Fixed

- Nucleic-acid inputs now reject incompatible scoring functions (`dfire`,
  `fastdfire`, `dfire2`, `pydock`, etc.) with a clear error and exit code 2.
  Only `dna` and `ddna` are accepted for nucleic-acid structures.
- Docking-worker failures are propagated to the process exit status instead of
  looking like a successful run.
- `auto` scoring can be resolved from the actual receptor/ligand residue names;
  an explicit incompatible selection is not silently replaced.

## Validation

- `cargo test --release --lib`: **29 passed**.
- Real protein–DNA regression (`MYB46-ALL.pdb` + `proCDH653-1137.pdb`,
  `score ... dna`): **-5194.963695**, not zero.
- A 1 swarm × 20 glowworms × 20 steps DNA docking produced a varying ranked
  score range (40.38–49.86) on Linux, confirming a non-constant scoring path.
  The same score was reproduced on macOS (`-5194.963695`).
- RNA/DNA PDB atom naming remains normalized by LKDock before docking.

## Binaries

- `LKlight-v1.2.1-linux-x86_64.tar.gz` — Linux x86-64
- `LKlight-v1.2.1-mac-arm64.zip` — macOS Apple Silicon
- `LKlight-v1.2.1-win-x64.zip` — Windows x64
