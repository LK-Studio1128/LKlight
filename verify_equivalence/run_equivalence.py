#!/usr/bin/env python3
"""
run_equivalence.py — LKlight vs Python LightDock accuracy-equivalence benchmark.

Pipeline per complex (Protein-Protein Docking Benchmark 5, lightdock_bm5):
  1. shared inputs  : `LKlight setup` generates setup.json + initial_positions_*.dat
  2. LKlight engine : `LKlight run` per swarm, then `LKlight generate` (all glowworms)
  3. Python engine  : `lightdock3` + `lgd_generate_conformations.py`
                     (or --py-mode official: parse official results/*.list as baseline)
  4. evaluation     : rank by scoring, CAPRI-class each top-N model vs native,
                      success rate = fraction of top-N with >=1 acceptable model

Outputs (under --outdir):
  summary.tsv             success rates per case & engine
  per_model_*.tsv         per-model CAPRI metrics
  success_curves.png      success-rate comparison plot (if matplotlib available)

Examples:
  python3 run_equivalence.py --bm5 /path/to/lightdock_bm5 --lklight target/release/LKlight
  python3 run_equivalence.py --fetch-bm5 ./bm5_data --n-cases 30 --swarms 25 --jobs 4
  python3 run_equivalence.py --bm5 ... --py-mode official     # official results as Python baseline
"""

import argparse
import concurrent.futures as cf
import csv
import glob
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

TOP_CUTOFFS = (1, 5, 10, 20, 50, 100)


def sh(cmd, cwd=None, timeout=7200, extra_env=None):
    """Run a command, returning (returncode, stdout)."""
    env = dict(os.environ)
    if extra_env:
        env.update(extra_env)
    try:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                           timeout=timeout, env=env)
        return r.returncode, r.stdout + r.stderr
    except Exception as exc:  # noqa: BLE001
        return -1, str(exc)


def find_tool(candidates, extra_paths=()):
    for name in candidates:
        for base in (sys.executable and os.path.dirname(sys.executable), *extra_paths):
            if base:
                p = os.path.join(base, name)
                if os.path.exists(p) and os.access(p, os.X_OK):
                    return p
        if shutil.which(name):
            return shutil.which(name)
    return None


def fetch_bm5(dest):
    if not os.path.isdir(os.path.join(dest, "data")):
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        rc, out = sh(["git", "clone", "--depth", "1",
                      "https://github.com/lightdock/lightdock_bm5.git", dest])
        if rc != 0:
            raise RuntimeError(f"git clone failed: {out[:500]}")
    return dest


def bm5_cases(bm5_root):
    data_dir = os.path.join(bm5_root, "data")
    if not os.path.isdir(data_dir):
        raise RuntimeError(f"no data/ dir in {bm5_root} — clone with --fetch-bm5")
    return sorted(d for d in os.listdir(data_dir)
                  if os.path.isdir(os.path.join(data_dir, d)))


def find_native(case_dir):
    for f in os.listdir(case_dir):
        if (f.endswith("_segid.pdb") or f.endswith(".segid.pdb")
                or f.endswith("_bound.pdb")):
            return os.path.join(case_dir, f)
    return None


# pdbtbx (used by LKlight) chokes on free-text REMARK lines such as
# "REMARK DATE:23-Dec-2018" (it tries to parse them as usize). Sanitize inputs
# by dropping metadata records; atom records are kept byte-for-byte.
_SANITIZE_KEEP = ("ATOM", "HETATM", "TER", "END", "CONECT", "CRYST1", "MODEL",
                  "MASTER", "HELIX", "SHEET", "SSBOND", "LINK", "SITE",
                  "ANISOU", "SEQRES", "TURN")


def sanitize_pdb(src, dst):
    kept = 0
    with open(src) as f, open(dst, "w") as g:
        for line in f:
            if line.startswith(_SANITIZE_KEEP):
                g.write(line)
                kept += 1
    if kept == 0:
        raise RuntimeError(f"sanitize_pdb: no atom records kept from {src}")
    return dst


def prepare_inputs(case_work, case_dir, lklight, swarms, glowworms, seed,
                   py_setup=None, py_pythonpath=None):
    """
    Prepare shared inputs with the OFFICIAL starting-point convention.

    LKlight's own setup places initial positions ~44 Å from the receptor centre,
    whereas the official LightDock setup samples surface points ~15-25 Å away;
    the latter is required for blind docking to converge. We therefore:
      1. sanitize the BM5 receptor/ligand PDBs,
      2. run `LKlight setup` once (produces the reference-frame structures
         lightdock_receptor_clean.pdb / lightdock_ligand_clean.pdb),
      3. run the official `lightdock3_setup.py` to get canonical
         init/initial_positions_*.dat (surface points),
      4. overwrite LKlight's own initial positions with the official ones and
         point setup.json's receptor_pdb/ligand_pdb at the plain clean names so
         the LKlight engine resolves `lightdock_<name>` correctly.
    """
    inputs = os.path.join(case_work, "inputs")
    os.makedirs(inputs, exist_ok=True)

    rec_candidates = glob.glob(os.path.join(case_dir, "*_A-noh.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*_A_noh.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*_rec*.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*receptor*.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*_A.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*.A.pdb"))
    lig_candidates = glob.glob(os.path.join(case_dir, "*_B-noh.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*_B_noh.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*_lig*.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*ligand*.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*_B.pdb")) + \
                     glob.glob(os.path.join(case_dir, "*.B.pdb"))
    if not rec_candidates or not lig_candidates:
        return None, None, None, inputs
    rec = rec_candidates[0]
    lig = lig_candidates[0]

    # 1) sanitize (pdbtbx chokes on free-text REMARK lines)
    rec_clean = os.path.join(inputs, "receptor_clean.pdb")
    lig_clean = os.path.join(inputs, "ligand_clean.pdb")
    sanitize_pdb(rec, rec_clean)
    sanitize_pdb(lig, lig_clean)

    # 2) LKlight setup -> reference-frame structures + a setup.json
    rc, out = sh([lklight, "setup", rec_clean, lig_clean, "-s", str(swarms),
                  "-g", str(glowworms), "--seed", str(seed)], cwd=inputs)
    if rc != 0:
        raise RuntimeError(f"LKlight setup failed: {out[:500]}")

    # 3) official starting points (surface points), if the official setup is usable
    if py_setup:
        official_dir = os.path.join(case_work, "inputs_official")
        os.makedirs(official_dir, exist_ok=True)
        # lightdock3_setup writes lightdock_*.pdb next to the INPUT PDBs, so give
        # it copies inside official_dir to avoid polluting the BM5 source data
        rec_copy = os.path.join(official_dir, "receptor_input.pdb")
        lig_copy = os.path.join(official_dir, "ligand_input.pdb")
        shutil.copy(rec, rec_copy)
        shutil.copy(lig, lig_copy)
        extra = {}
        shim = os.path.join(os.path.dirname(os.path.abspath(__file__)), "numpy_compat")
        extra["PYTHONPATH"] = os.pathsep.join(
            [p for p in (shim, py_pythonpath or "", os.environ.get("PYTHONPATH", "")) if p])
        cmd = [sys.executable, py_setup, rec_copy, lig_copy,
               "-s", str(swarms), "-g", str(glowworms),
               "--seed_points", str(seed), "--seed_anm", str(seed), "--noxt"]
        rc, out = sh(cmd, cwd=official_dir, extra_env=extra)
        if rc == 0 and os.path.isdir(os.path.join(official_dir, "init")):
            # 4) use the official surface points; keep LKlight's processed PDBs
            for s in range(swarms):
                src = os.path.join(official_dir, "init", f"initial_positions_{s}.dat")
                if os.path.exists(src):
                    shutil.copy(src, os.path.join(inputs, f"initial_positions_{s}.dat"))
            setup_path = os.path.join(inputs, "setup.json")
            cfg = common.read_setup(setup_path)
            cfg["receptor_pdb"] = "receptor_clean.pdb"
            cfg["ligand_pdb"] = "ligand_clean.pdb"
            with open(setup_path, "w") as fh:
                json.dump(cfg, fh, indent=2)
        else:
            print(f"  [prepare] official setup unavailable, using LKlight's own "
                  f"initial positions ({out[-200:]})", flush=True)

    return rec_clean, lig_clean, os.path.join(inputs, "setup.json"), inputs


def copy_inputs(work, inputs):
    """Copy setup.json + initial_positions + setup-generated lightdock_*.pdb into work."""
    for f in os.listdir(inputs):
        if f.startswith("lightdock_") and f.endswith(".pdb"):
            shutil.copy(os.path.join(inputs, f), work)
    for f in ("setup.json",):
        src = os.path.join(inputs, f)
        if os.path.exists(src):
            shutil.copy(src, work)


def setup_paths(setup_file):
    """Receptor/ligand paths as recorded in setup.json (relative to its dir)."""
    cfg = common.read_setup(setup_file)
    return cfg.get("receptor_pdb"), cfg.get("ligand_pdb")


def processed_pdb(work_dir, setup_file, which):
    """
    Resolve the PDB actually used by the engine inside work_dir: LightDock
    convention is `lightdock_<original_basename>.pdb` next to setup.json.
    """
    rec, lig = setup_paths(setup_file)
    base = os.path.basename(rec if which == "rec" else (lig or ""))
    if base:
        cand = os.path.join(work_dir, "lightdock_" + base)
        if os.path.exists(cand):
            return cand
    hits = sorted(glob.glob(os.path.join(work_dir, "lightdock_*.pdb")))
    return hits[0] if hits else None


def run_lklight(case_work, lklight, inputs, swarms, glowworms, steps, scoring):
    lk = os.path.join(case_work, "lk")
    os.makedirs(lk, exist_ok=True)
    copy_inputs(lk, inputs)
    for s in range(swarms):
        src = os.path.join(inputs, f"initial_positions_{s}.dat")
        if os.path.exists(src):
            shutil.copy(src, lk)
        os.makedirs(os.path.join(lk, f"swarm_{s}"), exist_ok=True)
    rec, lig = processed_pdb(lk, os.path.join(lk, "setup.json"), "rec"), \
               processed_pdb(lk, os.path.join(lk, "setup.json"), "lig")
    if not rec or not lig:
        return {"error": "processed receptor/ligand PDB missing for generate"}
    # run
    for s in range(swarms):
        pos = os.path.join(inputs, f"initial_positions_{s}.dat")
        if not os.path.exists(pos):
            continue
        rc, out = sh([lklight, "run", "setup.json", pos, str(steps), scoring], cwd=lk)
        if rc != 0:
            return {"error": f"lk run swarm {s}: {out[:300]}"}
    # generate all conformations per swarm
    for s in range(swarms):
        gso = os.path.join(lk, f"swarm_{s}", f"gso_{steps}.out")
        if not os.path.exists(gso):
            continue
        rc, out = sh([lklight, "generate", rec, lig, gso, str(glowworms + 1)], cwd=lk)
        if rc != 0:
            return {"error": f"lk generate swarm {s}: {out[:300]}"}
    return {}


def run_python(case_work, py_lightdock, py_generate, inputs, swarms,
               glowworms, steps, scoring, py_cores, py_pythonpath=None, py_data=None):
    py = os.path.join(case_work, "py")
    os.makedirs(py, exist_ok=True)
    copy_inputs(py, inputs)
    # Python LightDock reads initial positions from ./init/ (DEFAULT_POSITIONS_FOLDER)
    init_dir = os.path.join(py, "init")
    os.makedirs(init_dir, exist_ok=True)
    for s in range(swarms):
        src = os.path.join(inputs, f"initial_positions_{s}.dat")
        if os.path.exists(src):
            shutil.copy(src, init_dir)
        os.makedirs(os.path.join(py, f"swarm_{s}"), exist_ok=True)
    # NOTE: the Python generate script prepends the `lightdock_` prefix itself,
    # so pass the plain names from setup.json (opposite of the LKlight CLI).
    rec, lig = setup_paths(os.path.join(py, "setup.json"))
    if not rec or not lig:
        return {"error": "receptor/ligand names missing in py setup.json"}
    # the plain-named originals are also needed (valid_file checks the argument)
    for _name in (rec, lig):
        _src = os.path.join(inputs, _name)
        if os.path.exists(_src):
            shutil.copy(_src, py)

    extra = {}
    if py_pythonpath:
        cur = os.environ.get("PYTHONPATH", "")
        # numpy-compat shim first so sitecustomize patches numpy aliases at startup
        shim = os.path.join(os.path.dirname(os.path.abspath(__file__)), "numpy_compat")
        extra["PYTHONPATH"] = os.pathsep.join(
            [p for p in (shim, py_pythonpath, cur) if p])
    if py_data:
        extra["LIGHTDOCK_DATA"] = py_data

    # run through an explicit interpreter: the lightdock bin scripts carry a
    # stale shebang (python3.9) that does not exist in this environment
    runner = [sys.executable]

    cmd = runner + [py_lightdock, "setup.json", str(steps), "-s", scoring]
    if py_cores and py_cores > 1:
        cmd += ["-c", str(py_cores)]
    rc, out = sh(cmd, cwd=py, extra_env=extra)
    if rc != 0:
        return {"error": f"lightdock3 run: {out[:400]}"}

    # lightdock3.py swallows exceptions and may exit 0; verify real output
    missing = [s for s in range(swarms)
               if not os.path.exists(os.path.join(py, f"swarm_{s}", f"gso_{steps}.out"))]
    if missing:
        return {"error": f"lightdock3 produced no gso_{steps}.out for swarms {missing}; {out[-400:]}"}

    for s in range(swarms):
        gso = os.path.join(py, f"swarm_{s}", f"gso_{steps}.out")
        if not os.path.exists(gso):
            continue
        rc, out = sh(runner + [py_generate, rec, lig, gso, str(glowworms + 1)],
                     cwd=py, extra_env=extra)
        if rc != 0:
            return {"error": f"py generate swarm {s}: {out[:300]}"}
    return {}


def official_baseline(case_dir, bm5_root):
    """Locate and parse official results/*.list (BLIND scenario) as the Python baseline."""
    patterns = [
        os.path.join(bm5_root, "results", case_dir, "BLIND.list"),
        os.path.join(bm5_root, "results", f"{case_dir}_BLIND.list"),
        os.path.join(bm5_root, "results", f"{case_dir}-BLIND.list"),
        os.path.join(bm5_root, "results", case_dir, "*BLIND*.list"),
        os.path.join(bm5_root, "results", f"{case_dir}*BLIND*.list"),
    ]
    for pat in patterns:
        hits = glob.glob(pat)
        if hits:
            return common.parse_official_list(hits[0])
    return None


def process_case(case, args):
    case_dir = os.path.join(args.bm5, "data", case)
    out_dir = os.path.join(args.outdir, case)
    os.makedirs(out_dir, exist_ok=True)
    result = {"case": case, "lk": None, "py": None}

    native = find_native(case_dir)
    if not native:
        return {**result, "error": "no native structure found"}

    try:
        rec, lig, setup_path, inputs = prepare_inputs(
            out_dir, case_dir, args.lklight, args.swarms, args.glowworms, args.seed,
            getattr(args, "py_setup", None), getattr(args, "py_pythonpath", None))
        if not rec or not lig:
            return {**result, "error": "receptor/ligand not found"}

        lk_res = run_lklight(out_dir, args.lklight, inputs,
                             args.swarms, args.glowworms, args.steps, args.scoring)
        if lk_res.get("error"):
            return {**result, "error": lk_res["error"]}
        lk_dir = os.path.join(out_dir, "lk")
        rec_lk = processed_pdb(lk_dir, os.path.join(lk_dir, "setup.json"), "rec")
        lig_lk = processed_pdb(lk_dir, os.path.join(lk_dir, "setup.json"), "lig")
        if not rec_lk or not lig_lk:
            return {**result, "error": "processed receptor/ligand PDB missing in lk dir"}
        result["lk"] = common.evaluate_engine(
            lk_dir, rec_lk, lig_lk, native, args.steps, args.swarms)

        if args.py_mode == "run":
            if not args.py_lightdock or not args.py_generate:
                return {**result, "error": "py_mode=run requires --py-lightdock and --py-generate"}
            py_res = run_python(out_dir, args.py_lightdock, args.py_generate, inputs,
                                args.swarms, args.glowworms, args.steps,
                                args.scoring, args.py_cores, args.py_pythonpath,
                                args.py_data)
            if py_res.get("error"):
                return {**result, "error": py_res["error"]}
            py_dir = os.path.join(out_dir, "py")
            rec_py = processed_pdb(py_dir, os.path.join(py_dir, "setup.json"), "rec")
            lig_py = processed_pdb(py_dir, os.path.join(py_dir, "setup.json"), "lig")
            if not rec_py or not lig_py:
                return {**result, "error": "processed receptor/ligand PDB missing in py dir"}
            result["py"] = common.evaluate_engine(
                py_dir, rec_py, lig_py, native, args.steps, args.swarms)
        elif args.py_mode == "official":
            result["py"] = official_baseline(case, args.bm5)
            if result["py"] is None:
                return {**result, "error": "official BLIND.list not found"}
        else:  # skip
            result["py"] = None
    except Exception as exc:  # noqa: BLE001
        return {**result, "error": str(exc)}

    # persist per-model detail
    for eng, data in (("lk", result["lk"]), ("py", result["py"])):
        if not data:
            continue
        with open(os.path.join(out_dir, f"per_model_{eng}.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["rank", "swarm", "glowworm", "scoring", "fnat", "lrmsd", "rec_rmsd", "capri_class"])
            for m in data["models"]:
                w.writerow([m["rank"], m["swarm"], m["glowworm"],
                            f"{m['scoring']:.6f}", f"{m['fnat']:.4f}",
                            f"{m['lrmsd']:.4f}", f"{m['rec_rmsd']:.4f}",
                            m["capri_class"]])
    return result


def write_summary(results, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        header = ["case", "engine", "n_solutions"] + [f"top{n}" for n in TOP_CUTOFFS] + ["status"]
        w.writerow(header)
        for r in results:
            status = r.get("error", "ok")
            for eng in ("lk", "py"):
                if r.get(eng) is None:
                    continue
                sr = r[eng]["success_rates"]
                w.writerow([r["case"], eng, r[eng]["n_solutions"]]
                           + [f"{sr.get(n, 0.0):.4f}" for n in TOP_CUTOFFS] + [status])
            if r.get("error") and r.get("lk") is None and r.get("py") is None:
                w.writerow([r["case"], "-", 0] + ["-"] * len(TOP_CUTOFFS) + [status])


def plot_summary(results, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # noqa: BLE001
        return False
    xs = list(TOP_CUTOFFS)
    fig, ax = plt.subplots(figsize=(7, 5))
    for eng, color, label in (("lk", "#185FA5", "LKlight (this work)"),
                              ("py", "#0F6E56", "Python LightDock")):
        vals = []
        for n in xs:
            ok, tot = 0, 0
            for r in results:
                if r.get(eng) and not r.get("error"):
                    ok += int(r[eng]["success_rates"].get(n, 0.0) * n)
                    tot += n
            vals.append(ok / tot if tot else 0.0)
        ax.plot(xs, vals, marker="o", color=color, label=label)
    ax.set_xlabel("Top-N")
    ax.set_ylabel("Success rate (>=1 acceptable model)")
    ax.set_title("LKlight vs Python LightDock — success-rate equivalence")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return True


def main():
    ap = argparse.ArgumentParser(description="LKlight vs Python LightDock equivalence benchmark")
    ap.add_argument("--bm5", help="lightdock_bm5 repo path")
    ap.add_argument("--fetch-bm5", help="clone lightdock_bm5 into this dir first")
    ap.add_argument("--lklight", default="target/release/LKlight",
                    help="LKlight binary (default target/release/LKlight)")
    ap.add_argument("--py-lightdock", default=None, help="lightdock3 executable")
    ap.add_argument("--py-generate", default=None, help="lgd_generate_conformations.py executable")
    ap.add_argument("--py-pythonpath", default=None,
                    help="PYTHONPATH for the Python engine (e.g. a lightdock source checkout)")
    ap.add_argument("--py-setup", default=None,
                    help="lightdock3_setup.py path (official starting points; "
                         "defaults to <py-pythonpath>/bin/lightdock3_setup.py)")
    ap.add_argument("--py-data", default=None,
                    help="LIGHTDOCK_DATA dir for the Python engine (DFIRE params; "
                         "defaults to <repo>/data if present)")
    ap.add_argument("--py-mode", choices=["run", "official", "skip"], default="run",
                    help="run = execute Python engine; official = parse official .list; skip = LKlight only")
    ap.add_argument("--py-cores", type=int, default=1)
    ap.add_argument("--n-cases", type=int, default=20, help="number of BM5 cases (default 20)")
    ap.add_argument("--systems", nargs="*", default=[], help="explicit case list (overrides --n-cases)")
    ap.add_argument("--swarms", type=int, default=10, help="swarms per case (official: 400; keep both engines equal)")
    ap.add_argument("--glowworms", type=int, default=100, help="glowworms per swarm (official: 200)")
    ap.add_argument("--steps", type=int, default=100, help="optimization steps (official: 100)")
    ap.add_argument("--scoring", default="fastdfire", help="scoring method (official uses fastdfire)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="equivalence_results")
    ap.add_argument("--jobs", type=int, default=2, help="parallel cases")
    args = ap.parse_args()

    if args.fetch_bm5:
        args.bm5 = fetch_bm5(args.fetch_bm5)
    if not args.bm5 or not os.path.isdir(os.path.join(args.bm5, "data")):
        ap.error("need --bm5 (or --fetch-bm5); lightdock_bm5 data dir not found")

    lklight = args.lklight if os.path.isabs(args.lklight) else os.path.abspath(args.lklight)
    if not os.path.exists(lklight):
        ap.error(f"LKlight binary not found: {lklight} (build with `cargo build --release`)")

    if args.py_mode == "run":
        if not args.py_lightdock:
            args.py_lightdock = find_tool(["lightdock3", "lightdock"])
        if not args.py_generate:
            args.py_generate = find_tool(["lgd_generate_conformations.py",
                                          "lgd_generate_conformations"])
        if not args.py_lightdock or not args.py_generate:
            ap.error("py_mode=run needs lightdock3 and lgd_generate_conformations.py "
                     "on PATH; or use --py-mode official")

    # official starting points are beneficial for every mode (LKlight included)
    if not args.py_setup and args.py_pythonpath:
        cand = os.path.join(args.py_pythonpath, "bin", "lightdock3_setup.py")
        if os.path.exists(cand):
            args.py_setup = cand

    os.makedirs(args.outdir, exist_ok=True)
    cases = args.systems or bm5_cases(args.bm5)[: args.n_cases]
    print(f"[equivalence] {len(cases)} cases | swarms={args.swarms} "
          f"glowworms={args.glowworms} steps={args.steps} scoring={args.scoring} "
          f"py-mode={args.py_mode}", flush=True)

    results = []
    if args.jobs > 1:
        with cf.ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(process_case, c, args): c for c in cases}
            for fut in cf.as_completed(futures):
                r = fut.result()
                results.append(r)
                tag = r.get("error", "ok")
                print(f"  [{r['case']}] {tag}", flush=True)
    else:
        for c in cases:
            r = process_case(c, args)
            results.append(r)
            print(f"  [{r['case']}] {r.get('error', 'ok')}", flush=True)

    write_summary(results, os.path.join(args.outdir, "summary.tsv"))
    ok = sum(1 for r in results if not r.get("error"))
    print(f"\n[equivalence] done: {ok}/{len(results)} cases OK -> {args.outdir}/summary.tsv")

    png = os.path.join(args.outdir, "success_curves.png")
    if plot_summary(results, png):
        print(f"[equivalence] plot saved: {png}")
    else:
        print("[equivalence] matplotlib not available; skipping plot")


if __name__ == "__main__":
    main()
