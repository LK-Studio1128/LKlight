#!/usr/bin/env bash
# Build lightdock for macOS
# Run this script on a macOS machine with Rust installed.
# Install Rust: https://rustup.rs/
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== LKlight macOS build ==="
echo "Rust version: $(rustc --version)"
echo "Cargo version: $(cargo --version)"

# Detect architecture
ARCH=$(uname -m)
echo "Host arch: $ARCH"

# Option A: native build (fastest)
cargo build --release
BINARY="$SCRIPT_DIR/target/release/LKlight"

# Option B (optional): universal binary for both arm64 and x86_64
# Uncomment the block below to build a fat binary instead:
# -----------------------------------------------------------
# rustup target add aarch64-apple-darwin x86_64-apple-darwin 2>/dev/null || true
# cargo build --release --target aarch64-apple-darwin
# cargo build --release --target x86_64-apple-darwin
# lipo -create \
#   target/aarch64-apple-darwin/release/LKlight \
#   target/x86_64-apple-darwin/release/LKlight  \
#   -output target/release/LKlight-universal
# BINARY="$SCRIPT_DIR/target/release/LKlight-universal"
# -----------------------------------------------------------

echo ""
echo "=== Build successful ==="
echo "Binary: $BINARY"
echo "Size:   $(du -sh "$BINARY" | cut -f1)"
file "$BINARY" || true

# Package: binary only (all data files are embedded in the binary)
DIST="$SCRIPT_DIR/dist/LKlight-mac"
rm -rf "$DIST"
mkdir -p "$DIST"
cp "$BINARY" "$DIST/LKlight"

echo ""
echo "=== Package ready ==="
echo "Output: $DIST/"
echo "  LKlight   (self-contained binary, all data embedded)"
echo ""
echo "Subcommands:"
echo "  setup      <rec.pdb> <lig.pdb> <method> [-s N] [-g N] [--anm]"
echo "  run        <setup.json> <initial_positions.dat> <steps> <method>"
echo "  generate   <rec.pdb> <lig.pdb> <gso.out> <N>"
echo "  cluster    <gso.out>"
echo "  rank       <num_swarms> <steps>"
echo "  top        <ranking.list> <N>"
echo "  filter     <ranking.list> <restraints.list>"
echo "  gso_to_csv <ranking.list> <out.csv>"
echo "  score      <rec.pdb> <lig.pdb> <method> [--tx X --ty Y --tz Z]"
echo "  trajectory <rec.pdb> <lig.pdb> <swarm_id> <glowworm_id> <steps>"
echo "  map_contacts <rec.pdb> <lig.pdb> <gso.out>"
echo "  pipeline   <rec.pdb> <lig.pdb> <method> [--threads N]"
echo ""
echo "Methods: dfire fastdfire dfire2 dna mj3h pydock cpydock sd vdw pisa sipper tobi ddna"
echo ""
echo "Quick test:"
echo "  cd <lightdock_work_dir>"
echo "  $DIST/LKlight setup rec.pdb lig.pdb dfire -s 10 -g 50"
echo "  $DIST/LKlight run setup.json initial_positions_0.dat 100 dfire"
