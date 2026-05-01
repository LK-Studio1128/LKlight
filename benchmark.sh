#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# LightDock Performance Benchmark
# Compares: Python (lightdock3) | Rust-orig (lightdock-rust) | LKlight (this work)
#
# CLI interfaces:
#   Python    : lightdock3 setup.json <steps> -s <scoring> -l <swarm_id>
#   Rust-orig : lightdock-rust setup.json initial_positions_<N>.dat <steps> <scoring>
#   Rust-opt  : LKlight run setup.json initial_positions_<N>.dat <steps> <scoring>
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
ORIG_DIR="${LIGHTDOCK_REF_DIR:-}"
NEW_RUST="${LKLIGHT_BIN:-$REPO_DIR/target/release/LKlight}"
ORIG_RUST="$ORIG_DIR/lightdock-rust"
PYTHON3_BIN="$ORIG_DIR/lightdock3"

STEPS=100
REPEATS=3   # average over 3 runs

BENCH_TMP="$REPO_DIR/bench_tmp"
mkdir -p "$BENCH_TMP"

BASE1="${LKLIGHT_BENCH_1PPE:-$REPO_DIR/example/1ppe}"
BASE2="${LKLIGHT_BENCH_1AZP:-$REPO_DIR/example/1azp}"

if [[ ! -x "$NEW_RUST" ]]; then
    echo "ERROR: LKlight binary not found or not executable: $NEW_RUST"
    echo "Build first with: cargo build --release"
    exit 1
fi

if [[ -z "$ORIG_DIR" || ! -x "$ORIG_RUST" || ! -x "$PYTHON3_BIN" ]]; then
    echo "ERROR: reference LightDock binaries not found."
    echo "Set LIGHTDOCK_REF_DIR to a directory containing:"
    echo "  lightdock-rust"
    echo "  lightdock3"
    exit 1
fi

if [[ ! -d "$BASE1" || ! -d "$BASE2" ]]; then
    echo "ERROR: benchmark fixture directories not found."
    echo "Set these environment variables to prepared LightDock benchmark cases:"
    echo "  LKLIGHT_BENCH_1PPE=/path/to/example/1ppe"
    echo "  LKLIGHT_BENCH_1AZP=/path/to/example/1azp"
    exit 1
fi

# ─── timed run: returns avg wall-clock ms over REPEATS ───────────────────────
run_timed() {
    local workdir="$1"; shift
    local cmd=("$@")
    local total=0
    for _ in $(seq 1 $REPEATS); do
        local t0=$( date +%s%N )
        (cd "$workdir" && "${cmd[@]}" > /dev/null 2>&1) || true
        local t1=$( date +%s%N )
        total=$(( total + (t1 - t0) / 1000000 ))
    done
    echo $(( total / REPEATS ))
}

# ─── copy example to fresh tmp, strip pre-existing GSO outputs ───────────────
prep() {
    local src="$1" dst="$2"
    rm -rf "$dst"; cp -r "$src" "$dst"
    find "$dst" -name 'gso_*.out' -delete 2>/dev/null || true
}

# ─── speedup helper ──────────────────────────────────────────────────────────
spd() { awk "BEGIN{printf \"%.1f\", $1/$2}"; }

echo "======================================================="
echo " LightDock Benchmark  (steps=$STEPS, repeats=$REPEATS)"
echo " Single swarm (swarm_0) per invocation for fair compare"
echo "======================================================="
echo ""

run_case() {
    local name="$1" base="$2" scoring="$3" swarm_id="${4:-0}"
    local posdat="initial_positions_${swarm_id}.dat"

    echo "─── ${name}  [${scoring}] ───────────────────────────────────"

    # Python: patch swarms→1 so it doesn't require all N initial_positions files
    prep "$base" "$BENCH_TMP/${name}_py"
    # replace swarms count with 1 in setup.json
    sed -i.bak 's/"swarms": *[0-9]*/"swarms": 1/' "$BENCH_TMP/${name}_py/setup.json"
    T_PY=$( run_timed "$BENCH_TMP/${name}_py" \
        "$PYTHON3_BIN" setup.json $STEPS -s "$scoring" )
    echo "  Python      : ${T_PY} ms"

    # Rust-orig
    prep "$base" "$BENCH_TMP/${name}_orig"
    T_OR=$( run_timed "$BENCH_TMP/${name}_orig" \
        "$ORIG_RUST" setup.json "$posdat" $STEPS "$scoring" )
    echo "  Rust-orig   : ${T_OR} ms"

    # Rust-opt
    prep "$base" "$BENCH_TMP/${name}_opt"
    T_OPT=$( run_timed "$BENCH_TMP/${name}_opt" \
        "$NEW_RUST" run setup.json "$posdat" $STEPS "$scoring" )
    echo "  Rust-opt    : ${T_OPT} ms"

    local s_or=$( spd $T_PY  $T_OR  )
    local s_op=$( spd $T_PY  $T_OPT )
    local s_oo=$( spd $T_OR  $T_OPT )
    echo "  Speedup  Rust-orig  vs Python   : ${s_or}×"
    echo "  Speedup  Rust-opt   vs Python   : ${s_op}×"
    echo "  Speedup  Rust-opt   vs Rust-orig: ${s_oo}×"
    echo ""

    # store for summary
    eval "T_PY_${name}=$T_PY"
    eval "T_OR_${name}=$T_OR"
    eval "T_OPT_${name}=$T_OPT"
    eval "S_OP_${name}=$s_op"
    eval "S_OO_${name}=$s_oo"
}

# Case 1 – 1PPE pydock (protein-protein, no ANM, large system)
run_case "1PPE_pydock"  "$BASE1" "pydock"

# Case 2 – 1PPE dfire (DFIRE hash-table scoring)
run_case "1PPE_dfire"   "$BASE1" "dfire"

# Case 3 – 1AZP dna + ANM
run_case "1AZP_dna"     "$BASE2" "dna"

# Case 4 – 1PPE cpydock (complex scoring with desolvation)
run_case "1PPE_cpydock" "$BASE1" "cpydock"

# ─── Summary table ────────────────────────────────────────────────────────────
echo "======================================================="
echo " SUMMARY  (swarm_0, ${STEPS} steps, avg of ${REPEATS} runs)"
echo "======================================================="
printf "%-24s %9s %9s %9s %10s %11s\n" \
    "Case" "Py(ms)" "Orig(ms)" "Opt(ms)" "Opt/Py×" "Opt/Orig×"
printf "%-24s %9s %9s %9s %10s %11s\n" \
    "1PPE pydock"  "$T_PY_1PPE_pydock"  "$T_OR_1PPE_pydock"  "$T_OPT_1PPE_pydock"  "${S_OP_1PPE_pydock}×"  "${S_OO_1PPE_pydock}×"
printf "%-24s %9s %9s %9s %10s %11s\n" \
    "1PPE dfire"   "$T_PY_1PPE_dfire"   "$T_OR_1PPE_dfire"   "$T_OPT_1PPE_dfire"   "${S_OP_1PPE_dfire}×"   "${S_OO_1PPE_dfire}×"
printf "%-24s %9s %9s %9s %10s %11s\n" \
    "1AZP dna+ANM" "$T_PY_1AZP_dna"     "$T_OR_1AZP_dna"     "$T_OPT_1AZP_dna"     "${S_OP_1AZP_dna}×"     "${S_OO_1AZP_dna}×"
printf "%-24s %9s %9s %9s %10s %11s\n" \
    "1PPE cpydock" "$T_PY_1PPE_cpydock" "$T_OR_1PPE_cpydock" "$T_OPT_1PPE_cpydock" "${S_OP_1PPE_cpydock}×" "${S_OO_1PPE_cpydock}×"

# ─── JSON for paper ───────────────────────────────────────────────────────────
cat > "$BENCH_TMP/results.json" <<JSON
{
  "platform": "macOS arm64 (Apple Silicon)",
  "steps": $STEPS, "repeats": $REPEATS, "swarms": 1,
  "glowworms_per_swarm": 200,
  "1ppe_pydock":  { "py_ms": $T_PY_1PPE_pydock,  "orig_ms": $T_OR_1PPE_pydock,  "opt_ms": $T_OPT_1PPE_pydock,  "opt_vs_py": "$S_OP_1PPE_pydock",  "opt_vs_orig": "$S_OO_1PPE_pydock"  },
  "1ppe_dfire":   { "py_ms": $T_PY_1PPE_dfire,   "orig_ms": $T_OR_1PPE_dfire,   "opt_ms": $T_OPT_1PPE_dfire,   "opt_vs_py": "$S_OP_1PPE_dfire",   "opt_vs_orig": "$S_OO_1PPE_dfire"   },
  "1azp_dna_anm": { "py_ms": $T_PY_1AZP_dna,     "orig_ms": $T_OR_1AZP_dna,     "opt_ms": $T_OPT_1AZP_dna,     "opt_vs_py": "$S_OP_1AZP_dna",     "opt_vs_orig": "$S_OO_1AZP_dna"     },
  "1ppe_cpydock": { "py_ms": $T_PY_1PPE_cpydock, "orig_ms": $T_OR_1PPE_cpydock, "opt_ms": $T_OPT_1PPE_cpydock, "opt_vs_py": "$S_OP_1PPE_cpydock", "opt_vs_orig": "$S_OO_1PPE_cpydock" }
}
JSON
echo ""
echo "JSON: $BENCH_TMP/results.json"
echo "Done."
