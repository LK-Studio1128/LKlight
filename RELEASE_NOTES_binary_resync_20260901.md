# LKlight v1.1.0 binary resync

**Sync the prebuilt binaries with the v1.1.0 source (atom filtering enabled).**

The v1.1.0 source tree already contains `apply_setup_atom_filters`
(`--noh` / `--noxt` / `--now` handling in `setup`), but the binaries shipped in
`LKlight-{macos-arm64,win64,linux86}/` were built on 2026-06-01 and **do not**
contain that code — so `--now` (remove water), `--noxt` (remove OXT) and `--noh`
(remove hydrogens) were silently ignored in the field.

This commit rebuilds all three platforms from the current `main` source and
replaces the packaged binaries:

| Platform | Old (Jun 1) | New (rebuilt) |
|---|---|---|
| macOS arm64 | `LKlight-macos-arm64/LKlight` | aarch64, native |
| Windows x64 | `LKlight-win64/LKlight.exe` | PE32+ console |
| Linux x86-64 | `LKlight-linux86/LKlight` | musl static-PIE |

## Verified

- `strings <binary> | grep -c "Atom filters"` → 1 for all three rebuilt binaries
  (0 for the old Jun-1 ones).
- End-to-end on macOS: `setup rec.pdb lig.pdb dna -s 5 -g 3 --noxt --now` on a
  receptor containing 3 HOH + 1 OXT now reports
  `Atom filters [noh=false noxt=true now=true]: removed 4 receptor atoms` and
  the saved `lightdock_*.pdb` contains zero waters / zero OXT.

## Notes

- All data files are embedded in the binary — no external parameter files needed.
- Binaries are self-contained; `dist/` mirrors the same rebuilt artifacts.
