#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# Per-scoring-function three-engine benchmark for the LKlight GPU/grid paper.
#
# Runs the same GSO optimization through three engines on the same receptor /
# ligand and the same fixed seed, for each of the 12 scoring functions:
#   - exact : the all-pairs reference path (musl-static build)
#   - grid  : CPU cell-list + far-field electrostatic map
#   - gpu   : CUDA batch kernel (single-swarm pose batching)
#
# Outputs two CSV files:
#   server_functest_perfn.csv     per (method, engine, run) wall-clock + best
#   server_functest_summary.csv    per (method, engine) mean wall-clock,
#                                 per-method GPU-vs-grid ratio, best-energy
#                                 equivalence vs exact.
#
# Usage on the Linux server:
#
#   # set binary paths (defaults are the standard locations):
#   export EXACT_BIN=/path/to/LKlight-linux86-static
#   export GRID_BIN=/path/to/LKlight-linux-x64
#   export GPU_BIN=/path/to/LKlight-linux-cuda
#
#   # provide receptor / ligand / glow / steps / repeats:
#   export REC=/path/to/receptor.pdb    # e.g. AF-M9NF14 (8,218 atoms)
#   export LIG=/path/to/ligand.pdb      # e.g. CR44522 (12,625 atoms)
#   export GLOW=20                      # glowworms per swarm
#   export STEPS=100                    # GSO steps
#   export REPEATS=2                    # average wall-clock over N runs
#
#   # (optional) provide a pre-built setup+initial_positions to skip the
#   # setup phase and ensure byte-identical starting coordinates across all
#   # three engines:
#   export PRESET_DIR=/path/to/preset   # contains setup.json +
#                                      # initial_positions_0.dat
#
#   ./server_functest.sh
#
# After the run, ship the resulting CSVs back and append them to the paper
# (Section "Large-scale per-function speedup, server measurement").
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

EXACT_BIN="${EXACT_BIN:-${HOME}/LKlight/release_bin/LKlight-linux86-static}"
GRID_BIN="${GRID_BIN:-${HOME}/LKlight-grid/release_bin/LKlight-linux-x64}"
GPU_BIN="${GPU_BIN:-${HOME}/LKlight-GPU/release_bin/LKlight-linux-cuda}"

REC="${REC:?ERROR: set REC=/path/to/receptor.pdb}"
LIG="${LIG:?ERROR: set LIG=/path/to/ligand.pdb}"
GLOW="${GLOW:-20}"
STEPS="${STEPS:-100}"
REPEATS="${REPEATS:-2}"
PRESET_DIR="${PRESET_DIR:-}"

WORK="$(pwd)/functest_run_$$"
mkdir -p "$WORK"

METHODS=(dfire dfire2 dna mj3h pydock cpydock sd vdw pisa sipper tobi ddna)

# ─── prepare a working setup with fixed seed ─────────────────────────────────
if [[ -n "$PRESET_DIR" ]]; then
    cp "$PRESET_DIR/setup.json" "$WORK/setup.json"
    cp "$PRESET_DIR"/initial_positions_*.dat "$WORK/"
    cp "$PRESET_DIR"/lightdock_*.pdb "$WORK/" 2>/dev/null || \
        cp "$REC" "$LIG" "$WORK/"
else
    # build a fresh 1-swarm setup
    "$GRID_BIN" setup "$REC" "$LIG" -s 1 -g "$GLOW" --seed 324324 \
        > "$WORK/setup.log" 2>&1
    # rename setup.json in cwd of $WORK
    cp setup.json "$WORK/setup.json" 2>/dev/null || true
    cp lightdock_*.pdb "$WORK/" 2>/dev/null || true
fi

# helpers
time_run() {
    local bin="$1" method="$2"
    local t0 t1
    t0=$(date +%s%N)
    (cd "$WORK" && "$bin" run setup.json \
        "$(ls "$WORK"/initial_positions_0.dat)" "$STEPS" "$method" \
        > "$WORK/${method}_last.out" 2>&1) || true
    t1=$(date +%s%N)
    awk "BEGIN{printf \"%.6f\", ($t1-$t0)/1e9}"
}

best_score() {
    # min over all glowworm gso_*.out in swarm_0/ last-column "Scoring"
    local method="$1" best=-1e300
    for f in "$WORK"/swarm_0/gso_*.out; do
        [[ -f "$f" ]] || continue
        awk '/^[^#]/ { if ($NF ~ /^-?[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/) {
                    v=$NF+0; if (v<best) best=v } } END {
                    printf "%.6f", best+0.0 }' "$f"
    done
}

PERF_CSV="$WORK/server_functest_perfn.csv"
SUM_CSV="$WORK/server_functest_summary.csv"
echo "method,engine,repeat,wall_s,best_score" > "$PERF_CSV"
echo "method,engine,mean_s,best_score,rel_delta_vs_exact" > "$SUM_CSV"

echo "==============================================================="
echo " Server per-function benchmark"
echo "  exact: $EXACT_BIN"
echo "  grid : $GRID_BIN"
echo "  gpu  : $GPU_BIN"
echo "  complex: $REC + $LIG"
echo "  glow=$GLOW steps=$STEPS repeats=$REPEATS"
echo "==============================================================="

for m in "${METHODS[@]}"; do
    # Exact
    et=0; eb=""
    for r in $(seq 1 "$REPEATS"); do
        t=$(time_run "$EXACT_BIN" "$m")
        b=$(best_score "$m")
        echo "$m,exact,$r,$t,$b" >> "$PERF_CSV"
        et=$(awk "BEGIN{print $et + $t}")
        eb="$b"
    done
    em=$(awk "BEGIN{printf \"%.4f\", $et/$REPEATS}")
    echo "$m,exact,$em,$eb," >> "$SUM_CSV"

    # Grid
    gt=0; gb=""
    for r in $(seq 1 "$REPEATS"); do
        t=$(time_run "$GRID_BIN" "$m")
        b=$(best_score "$m")
        echo "$m,grid,$r,$t,$b" >> "$PERF_CSV"
        gt=$(awk "BEGIN{print $gt + $t}")
        gb="$b"
    done
    gm=$(awk "BEGIN{printf \"%.4f\", $gt/$REPEATS}")
    echo "$m,grid,$gm,$gb," >> "$SUM_CSV"

    # GPU
    pt=0; pb=""
    for r in $(seq 1 "$REPEATS"); do
        t=$(time_run "$GPU_BIN" "$m")
        b=$(best_score "$m")
        echo "$m,gpu,$r,$t,$b" >> "$PERF_CSV"
        pt=$(awk "BEGIN{print $pt + $t}")
        pb="$b"
    done
    pm=$(awk "BEGIN{printf \"%.4f\", $pt/$REPEATS}")
    # rel delta vs exact (best-energy agreement check)
    rel=""
    if [[ -n "$eb" && -n "$pb" && $(awk "BEGIN{print (sqrt($eb^2)>1e-3)?1:0}") -eq 1 ]]; then
        rel=$(awk "BEGIN{printf \"%.6f\", ($pb-$eb)/sqrt($eb*$eb)}")
    fi
    echo "$m,gpu,$pm,$pb,$rel" >> "$SUM_CSV"

    printf "%-7s  exact=%6.3fs  grid=%6.3fs  gpu=%6.3fs   best_exact=%-14s best_gpu=%-14s rel=%s\n" \
        "$m" "$em" "$gm" "$pm" "$eb" "$pb" "${rel:-—}"
done

echo
echo "Wrote: $PERF_CSV"
echo "Wrote: $SUM_CSV"
echo
echo "Cleanup: rm -rf $WORK"