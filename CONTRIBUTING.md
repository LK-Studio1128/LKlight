# Contributing to LKlight

Thank you for your interest in contributing to LKlight!

## License Agreement

LKlight is distributed under the **GNU General Public License v3.0 or later
(GPL-3.0-or-later)**. By submitting a Pull Request or any other contribution,
you agree that your work will be incorporated into the project and distributed
under the same license. You confirm that you have the right to make the
contribution under these terms.

LKlight is a derivative of LightDock (GPL-3.0). Please do not submit code
that was taken from sources with incompatible licenses (e.g., proprietary,
AGPL-only, or CC-BY-NC).

---

## Reporting Bugs

Open a [GitHub Issue](../../issues) and include:

1. **LKlight version** — output of `LKlight --version` or `git rev-parse --short HEAD`
2. **Operating system and architecture** (e.g., `macOS 14.4 arm64`, `Ubuntu 22.04 x86_64`)
3. **Rust toolchain** — output of `rustc --version`
4. **Minimal reproducible example** — PDB files + exact command line
5. **Expected output** vs. **actual output / error message**
6. If the issue involves numerical results, include the Python LightDock
   reference output for comparison

---

## Contributing Code

### 1. Set Up Your Development Environment

```bash
# Install Rust stable (rustup recommended)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Clone your fork
git clone https://github.com/<your-username>/LKlight.git
cd LKlight

# Verify the build
cargo build --release

# Run the test suite
cargo test --lib      # 29 unit tests must all pass
```

### 2. Create a Feature Branch

```bash
git checkout -b feature/my-improvement
# or
git checkout -b fix/dfire-edge-case
```

### 3. Coding Style

- Follow standard Rust idioms. Run `cargo fmt` before committing.
- Run `cargo clippy -- -D warnings` and resolve all warnings.
- Prefer **safe Rust**. Use `unsafe` only with a documented safety invariant comment.
- For fixed-size arrays known at compile time, prefer `[T; N]` (stack) over `Vec<T>` (heap).
- Public API items (structs, traits, public functions) should have `///` doc-comments.
- Do not add `println!` debug output; use the `log` crate (`log::debug!`, `log::info!`).

### 4. Adding or Modifying Scoring Functions

If you add or significantly modify a scoring function:

1. Add a corresponding unit test in the same file with a known-good reference
   value from Python LightDock (tolerance ε ≤ 10⁻⁸).
2. Verify the test passes: `cargo test --lib <function_name>`.
3. Document the scoring function's citation (journal paper) in a `///` comment
   at the top of the module.
4. If the modification changes numerical output by more than ε = 10⁻⁸,
   explain why in the PR description (e.g., improved numerical stability,
   corrected algorithm).

### 5. Performance Optimizations

For any optimization that changes runtime behavior:

1. Provide before/after benchmark numbers. Use `bash benchmark.sh` or a
   targeted `cargo bench` if a criterion benchmark exists.
2. Confirm all 29 unit tests and 160 integration tests still pass.
3. Explain the algorithmic rationale (why the optimization is correct) and
   quantify the expected gain.
4. Be aware of the **spatial grid lesson** (see `PAPER.md §4.3`): grids with
   cutoffs ≥ 15 Å typically regress performance for small protein systems.
   Always measure before and after.

### 6. Submit a Pull Request

1. Push your branch and open a Pull Request against `main`.
2. Describe:
   - The problem addressed or feature added
   - The approach taken
   - Benchmark data (if performance-related)
   - Any tradeoffs or limitations
3. CI must pass (build + `cargo test`).

---

## Code Review Process

- All PRs require at least one review before merging.
- Maintainers may request changes for correctness, style, or performance.
- Numerical changes to scoring functions will be held to a higher bar and
  require comparison against Python LightDock reference values.

---

## Attribution

All contributors are listed in `CHANGELOG.md`. The original LightDock
authors (Brian Jiménez-García et al.) are permanently attributed in `NOTICE`.

---

## Questions

Open an issue tagged **`question`** for general questions about the codebase
or algorithm. For LightDock algorithm questions, the [LightDock documentation](https://lightdock.org)
and original papers [1,2] are the primary references.

---

*[1] Jiménez-García et al., Bioinformatics 2018 — doi:10.1093/bioinformatics/btx555*  
*[2] Roel-Touris et al., Bioinformatics 2020 — doi:10.1093/bioinformatics/btz642*
