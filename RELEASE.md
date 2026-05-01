# LKlight Release Guide

This document describes how to publish LKlight source code and pre-built binaries.

## 1. Repository contents

The GitHub repository should include:

- `src/` — Rust source code
- `tests/` — lightweight test PDB fixtures
- `data/` — parameter files required for embedded scoring data
- `Cargo.toml` and `Cargo.lock` — reproducible Rust binary build
- `README.md` — user-facing overview and usage
- `LICENSE` — GPL-3.0-or-later license text
- `NOTICE` — upstream LightDock attribution and modification summary
- `CHANGELOG.md` — differences from upstream `lightdock-rust`
- `PAPER.md` — technical manuscript / method documentation
- `CONTRIBUTING.md` — contribution and license rules
- `.github/workflows/rust.yml` — CI build and test workflow
- `build_mac.sh`, `build_linux.sh`, `build_win.bat` — platform build helpers

The repository should not include generated build artifacts:

- `target/`
- `dist/`
- `bench_tmp/`
- `.DS_Store`
- backup files such as `*.bak`

## 2. Build checks before release

```bash
cargo fmt --check
cargo test --lib
cargo build --release
```

Optional lint:

```bash
cargo clippy -- -D warnings
```

## 3. Platform binary packages

Recommended release asset names:

| Platform | Asset name | Build command |
|---|---|---|
| macOS arm64 | `LKlight-macos-arm64.tar.gz` | `bash build_mac.sh` |
| Linux x86_64 | `LKlight-linux-x86_64.tar.gz` | `bash build_linux.sh` |
| Windows x64 | `LKlight-windows-x64.zip` | `build_win.bat` |

## 4. GitHub release notes template

```markdown
# LKlight v1.0.0

High-performance Rust implementation of the LightDock docking engine.

## Highlights

- 12 scoring-function families with 13 command-line method names: dfire, fastdfire, dfire2, dna, mj3h, pydock, cpydock, sd, vdw, pisa, sipper, tobi, ddna
- Single self-contained binary
- Fixes critical upstream `lightdock-rust` issues in DFIRE parameters, ANM stride, non-standard residues, and ANM atom-count assertion
- Rayon parallelization and SIMD-friendly hot loops
- macOS / Linux / Windows build helpers

## Benchmark summary

macOS arm64, swarm_0, 200 glowworms, 100 steps, 3-run average:

| Case | Python LightDock | lightdock-rust baseline | LKlight |
|---|---:|---:|---:|
| 1PPE pydock | 858 ms | 7693 ms | 290 ms |
| 1PPE dfire | 840 ms | crash | 33 ms |
| 1AZP dna+ANM | 760 ms | 14142 ms | 46 ms |
| 1PPE cpydock | 844 ms | 7158 ms | 44 ms |

## License

LKlight is a GPL-3.0-or-later derivative work of LightDock / lightdock-rust.
Please keep `LICENSE` and `NOTICE` when redistributing.
```

## 5. Compliance notes

LKlight is a derivative of LightDock and `lightdock-rust`, which are GPL-3.0-or-later projects. When publishing or redistributing LKlight:

- keep the GPL license text;
- keep upstream attribution in `NOTICE`;
- publish the modified source code when distributing binaries;
- do not mix LKlight code into proprietary closed-source products unless GPL obligations are satisfied.
