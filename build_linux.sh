#!/usr/bin/env bash
# Build lightdock-rust for Linux (native)
# Run this script on a Linux machine with Rust installed.
# Install Rust: curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== LKlight Linux build ==="
echo "Rust version: $(rustc --version)"
echo "Cargo version: $(cargo --version)"

# Ensure musl target is available for maximum portability (static binary)
# Comment out the next two lines if you prefer glibc dynamic linking
TARGET="x86_64-unknown-linux-musl"
rustup target add "$TARGET" 2>/dev/null || true

# Try musl static build first (most portable), fall back to native glibc
if rustup target list --installed 2>/dev/null | grep -q "x86_64-unknown-linux-musl"; then
    echo "Building static musl binary (portable across Linux distros)..."
    cargo build --release --target "$TARGET"
    BINARY="$SCRIPT_DIR/target/$TARGET/release/LKlight"
else
    echo "Building native glibc binary..."
    cargo build --release
    BINARY="$SCRIPT_DIR/target/release/LKlight"
fi

echo ""
echo "=== Build successful ==="
echo "Binary: $BINARY"
echo "Size:   $(du -sh "$BINARY" | cut -f1)"

# Verify it's actually a Linux ELF
file "$BINARY" || true

# Package: binary only (all data files are embedded in the binary)
DIST="$SCRIPT_DIR/dist/LKlight-linux"
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
