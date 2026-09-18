#!/bin/bash
# mac_drive_ht1000.sh — Mac mini M4 driver for the Table-4 1AZP 1000x1000 benchmark.
# Mirrors win_drive_gpu.ps1 (same seed 324324, same PDBs, same initial positions,
# same setup.json) so the per-engine results are directly comparable with the
# RTX 5090 same-batch run (functest/ht1000_summary_5090.csv).
#
# Usage: mac_drive_ht1000.sh <engine>    engine ∈ {exact, grid}
#   EXE  : absolute path to the v1.2.0 mac-arm64 binary (set per engine below).
#   WDIR : the engine's work tree (PDBs + setup.json + initial_positions_0.dat).
# Run engines SEQUENTIALLY (each saturates all 10 cores).
set -u
ENGINE="$1"
ROOT="$(cd "$(dirname "$0")" && pwd)"
case "$ENGINE" in
  exact) EXE="$ROOT/../release_v120/exact_out/LKlight-v1.2.0-mac-arm64/LKlight" ;;
  grid)  EXE="$ROOT/../release_v120/grid/LKlight-mac-arm64" ;;
  *) echo "usage: mac_drive_ht1000.sh exact|grid"; exit 2 ;;
esac
WDIR="$ROOT/$ENGINE"
[ -x "$EXE" ] || { echo "missing engine: $EXE"; exit 3; }

cd "$WDIR" || exit 4
rm -rf swarm_0

best=$(python3 - <<'PY'
import glob,sys
best=None
for f in glob.glob("swarm_0/gso_*.out"):
    try:
        for ln in open(f):
            s=ln.strip()
            if not s or s.startswith("#"): continue
            try: v=float(s.split()[-1])
            except ValueError: continue
            if best is None or v<best: best=v
    except FileNotFoundError: pass
print("" if best is None else repr(best))
PY
)

echo ">>> $ENGINE (1000 glow x 1000 steps, dna, 1azp) on $(sysctl -n machdep.cpu.brand_string) ..."
START=$(python3 -c 'import time;print(time.time())')
/usr/bin/time -p "$EXE" run setup.json initial_positions_0.dat 1000 dna > stdout.log 2> stderr.log
RC=$?
END=$(python3 -c 'import time;print(time.time())')
WALL=$(python3 -c "print(round($END-$START,2))")
# best scan (after run)
BEST=$(python3 - <<'PY'
import glob
best=None
for f in glob.glob("swarm_0/gso_*.out"):
    for ln in open(f):
        s=ln.strip()
        if not s or s.startswith("#"): continue
        try: v=float(s.split()[-1])
        except ValueError: continue
        if best is None or v<best: best=v
print("NA" if best is None else ("%.10f"%best))
PY
)
echo "rc=$RC wall=${WALL}s best=$BEST"
echo "$ENGINE,$RC,$WALL,$BEST,$(date '+%Y-%m-%d %H:%M')" >> "$ROOT/ht1000_summary_mac.csv"
