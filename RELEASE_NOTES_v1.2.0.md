# LKlight v1.2.0

**DNA scoring speedup + heavy-atom clash fix + auditing tools.** This release
synchronises the original all-pairs reference engine with the grid / GPU
accelerated engines (LKlight-grid v1.2.0, LKlight-GPU v1.2.0) and ships
rebuilt binaries for Linux x86_64, macOS arm64 and Windows x64.

## Changed

- **1-D Z-window DNA scoring speedup** (`src/dna*.rs`): the all-atom pairwise
  scan for the `dna` / `ddna` families prunes receptor atoms with a 1-D
  Z-window along the long receptor axis, reducing the enumerated pair set on
  long nucleic-acid ligands. Linux server rebuild (static-pie) verified as part
  of this release.

## Fixed

- **Systematic interpenetration in DNA scoring — heavy-atom clash penalty**
  (`src/dna*.rs`): poses burying ligand atoms into receptor backbone heavy
  atoms were not penalised, allowing deep clashes; a clash penalty term now
  corrects the interpenetration found in the clash-fix campaign (2026-09-02
  binary resync).

## Added

- **Clash / contact auditing subcommands** (e.g. `map_contacts`) for docking
  poses; multi-scenario acceptance fixtures and a Zenodo DOI badge;
  `--noh / --noxt / --now` atom filters enabled in prebuilt binaries.

## Binaries (v1.2.0)

- `LKlight-v1.2.0-linux-x86_64.tar.gz` — Ubuntu server build, static-pie (musl)
- `LKlight-v1.2.0-mac-arm64.zip` — macOS Apple Silicon
- `LKlight-v1.2.0-win-x64.zip` — Windows x64 (UCRT)

See `CHANGELOG.md` and `README.md` for details.
