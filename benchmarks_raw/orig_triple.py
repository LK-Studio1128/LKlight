#!/usr/bin/env python3
"""orig_triple.py — 31-decoy, three-way audit: original LightDock 0.9.4 (python)
vs LKlight exact vs LKlight-grid (v1.2.0, 0.5 A default).

Question being answered: is the "<=0.5%/bit-identical" claim *against the
original LightDock*, not merely against our own exact engine?

Design (file-level poses, copied from numeric_validation_v2.py): the ligand
PDB is rigidly transformed on disk first, then every scorer runs at identity.
All three sides therefore see byte-identical input; any diff is a real one.

Poses: native + 30 deterministic rigid-body decoys (same RNG/seed as the
§3.7 equiv_harness: translations 0.5-26 A, rotations <= 60 deg).

Outputs (into this dir): orig_triple_perpose.csv, orig_triple_summary.csv
"""
import sys, os, math, random, statistics, subprocess, csv, importlib, time

# ---- env must be pinned before any numpy/scipy import ----------------------
PY39_SP = "/Users/luoxiaowen/mambaforge/envs/docking_env/lib/python3.9/site-packages"
sys.path.insert(1, PY39_SP)  # numpy 1.23.5 + scipy 1.13.1 + lightdock 0.9.4 (cp39)
import numpy as np
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/Users/luoxiaowen/Desktop/LKDock/byi"
EXACT = f"{ROOT}/LKlight/target/release/LKlight"          # all-pairs exact
GRID  = f"{ROOT}/LKlight-grid/target/release/LKlight"      # CPU grid, v1.2.0 0.5A
WORK  = "/tmp/numval_orig"

SYSTEMS = {  # (receptor, ligand) — both already hydrogen-stripped (noh)
    "2oob": (f"{WORK}/2oob_receptor.pdb", f"{WORK}/2oob_ligand.pdb"),
    "1azp": (f"{WORK}/1azp_receptor_noh.pdb", f"{WORK}/1azp_ligand_noh.pdb"),
}
PROTEIN_METHODS = ["dfire", "fastdfire", "dfire2", "mj3h", "cpydock",
                   "sd", "vdw", "pisa", "sipper", "tobi"]
DNA_METHODS = ["dna", "ddna"]
# pydock has NO python-side module in LightDock 0.9.4 (pyDock family == cpydock);
# kept for a rust-only exact-vs-grid echo of the §3.7 audit.
PYDOCK_ONLY = ["pydock"]
NDECOYS = 30  # -> 31 poses incl. native


def normq(qw, qx, qy, qz):
    n = math.sqrt(qw*qw + qx*qx + qy*qy + qz*qz) or 1.0
    return qw/n, qx/n, qy/n, qz/n


def quat_from_axis(axis, deg):
    r = math.radians(deg); s = math.sin(r/2); c = math.cos(r/2)
    ax, ay, az = axis
    return normq(c, ax*s, ay*s, az*s)


def make_decoys(n):
    """Same generator/seed as equiv_harness.make_decoys (20260904)."""
    decoys = [(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)]
    rng = random.Random(20260904)
    for _ in range(n):
        r = rng.uniform(0, 1)
        mag = rng.uniform(2.0, 8.0) if r < 0.6 else rng.uniform(10.0, 26.0)
        th = rng.uniform(0, 2*math.pi); ph = math.acos(rng.uniform(-1, 1))
        tx = mag*math.sin(ph)*math.cos(th)
        ty = mag*math.sin(ph)*math.sin(th)
        tz = mag*math.cos(ph)
        deg = rng.uniform(0, 60)
        ax = [rng.uniform(-1, 1) for _ in range(3)]
        L = math.sqrt(sum(a*a for a in ax)) or 1.0
        ax = [a/L for a in ax]
        qw, qx, qy, qz = quat_from_axis(ax, deg)
        decoys.append((tx, ty, tz, qw, qx, qy, qz))
    return decoys


def quat_to_R(qw, qx, qy, qz):
    qw, qx, qy, qz = normq(qw, qx, qy, qz)
    return np.array([
        [1-2*(qy*qy+qz*qz), 2*(qx*qy-qw*qz),   2*(qx*qz+qw*qy)],
        [2*(qx*qy+qw*qz),   1-2*(qx*qx+qz*qz), 2*(qy*qz-qw*qx)],
        [2*(qx*qz-qw*qy),   2*(qy*qz+qw*qx),   1-2*(qx*qx+qy*qy)],
    ])


def write_transformed(src, dst, R=None, t=None):
    out = []
    for line in open(src):
        if line.startswith(("ATOM", "HETATM")):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if R is not None:
                xyz = R @ xyz
            if t is not None:
                xyz = xyz + t
            out.append(line[:30] + "%8.3f%8.3f%8.3f" % tuple(xyz) + line[54:])
        else:
            out.append(line)
    open(dst, "w").writelines(out)


# ---------------- python original-LightDock side -----------------------------
_cache = {}
def get_parsed(path):
    if "pose_" in os.path.basename(path):   # pose files: never cache (stale-risk)
        atoms, residues, chains = parse_complex_from_file(path)
        return Complex(chains, atoms, structure_file_name=path)
    if path not in _cache:
        atoms, residues, chains = parse_complex_from_file(path)
        _cache[path] = Complex(chains, atoms, structure_file_name=path)
    return _cache[path]

def py_score(method, rec_path, lig_path):
    receptor = get_parsed(rec_path).clone()
    ligand = get_parsed(lig_path).clone()
    module = importlib.import_module("lightdock.scoring.%s.driver" % method)
    ScoringFunction = getattr(module, "DefinedScoringFunction")
    ModelAdapter = getattr(module, "DefinedModelAdapter")
    adapter = ModelAdapter(receptor, ligand)
    sf = ScoringFunction()
    return float(sf(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                    adapter.ligand_model, adapter.ligand_model.coordinates[0]))


def coords(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path) if l.startswith(("ATOM", "HETATM"))])

def min_inter_dist(rec_c, lig_c):
    if len(rec_c) == 0 or len(lig_c) == 0:
        return float("nan")
    m = 1e30
    for i in range(0, len(rec_c), 256):
        d = ((rec_c[i:i+256, None, :] - lig_c[None, :, :])**2).sum(-1)
        m = min(m, d.min())
    return math.sqrt(m)


def rust_score(binpath, method, rec_path, lig_path):
    r = subprocess.run([binpath, "score", rec_path, lig_path, method,
                        "--tx", "0", "--ty", "0", "--tz", "0",
                        "--qw", "1", "--qx", "0", "--qy", "0", "--qz", "0"],
                       capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        return None, "rc=%d: %s" % (r.returncode, r.stderr.strip()[:120])
    for tok in reversed(r.stdout.replace(":", " ").replace(",", " ").split()):
        try:
            return float(tok), ""
        except ValueError:
            continue
    return None, "unparseable stdout"


def spearman(x, y):
    def rk(s):
        idx = sorted(range(len(s)), key=lambda i: s[i])
        r = [0]*len(s); i = 0
        while i < len(s):
            j = i
            while j+1 < len(s) and s[idx[j+1]] == s[idx[i]]:
                j += 1
            a = (i+j)/2
            for k in range(i, j+1):
                r[idx[k]] = a
            i = j+1
        return r
    n = len(x); mx = sum(x)/n; my = sum(y)/n
    rx, ry = rk(x), rk(y)
    cov = sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    vx = math.sqrt(sum((v-mx)**2 for v in rx)); vy = math.sqrt(sum((v-my)**2 for v in ry))
    return cov/(vx*vy) if vx*vy else float("nan")


def main():
    decoys = make_decoys(NDECOYS)
    os.makedirs(HERE, exist_ok=True)
    perpose = open(os.path.join(HERE, "orig_triple_perpose.csv"), "w", newline="")
    pw = csv.writer(perpose)
    pw.writerow(["pose", "method", "system", "python", "exact", "grid",
                 "diff_exact_py", "diff_grid_py", "rel_grid_py", "min_dist"])

    A = {}
    for m in PROTEIN_METHODS: A[m] = "2oob"
    for m in DNA_METHODS:     A[m] = "1azp"
    for m in PYDOCK_ONLY:     A[m] = "1azp"
    REC_COORDS = {sysname: coords(SYSTEMS[sysname][0]) for sysname in set(A.values())}
    LIG_COORDS = {sysname: coords(SYSTEMS[sysname][1]) for sysname in set(A.values())}

    t_start = time.time()
    n_ok = 0; n_items = 0
    for pidx, (tx, ty, tz, qw, qx, qy, qz) in enumerate(decoys):
        R = quat_to_R(qw, qx, qy, qz)
        t = np.array([tx, ty, tz])
        tag = "native" if pidx == 0 else "decoy%02d" % pidx
        for method, sysname in A.items():
            recf, ligf = SYSTEMS[sysname]
            pose_lig = os.path.join(WORK, "pose_%s_%s_ligand.pdb" % (tag, sysname))
            write_transformed(ligf, pose_lig, R=R, t=t)
            lig_c = (R @ LIG_COORDS[sysname].T).T + t
            md = min_inter_dist(REC_COORDS[sysname], lig_c)
            # python original
            try:
                ps = py_score(method, recf, pose_lig)
                perr = ""
            except Exception as e:
                ps = None
                perr = "%s: %s" % (type(e).__name__, str(e)[:100])
            # rust exact + grid (identity, same file)
            es, eerr = rust_score(EXACT, method, recf, pose_lig)
            gs, gerr = rust_score(GRID, method, recf, pose_lig)
            if ps is not None and es is not None and gs is not None:
                de = abs(es - ps); dg = abs(gs - ps)
                relg = dg/abs(ps) if abs(ps) > 1e-3 else dg
                if de < 1e-5:
                    n_ok += 1
                n_items += 1
                pw.writerow([tag, method, sysname, "%.10f" % ps,
                             "%.6f" % es, "%.6f" % gs,
                             "%.3g" % de, "%.3g" % dg, "%.4g" % relg,
                             "%.2f" % md])
                print("%-8s %-10s exact-py=%.2e grid-py=%.2e (rel %.3g, md %.1fA)" %
                      (tag, method, de, dg, relg, md), flush=True)
            else:
                pw.writerow([tag, method, sysname,
                             ps if ps is not None else "ERR:%s" % perr,
                             es if es is not None else "ERR:%s" % eerr,
                             gs if gs is not None else "ERR:%s" % gerr,
                             "", "", "", "%.2f" % md])
                print("%-8s %-10s ERROR py=%s exact=%s grid=%s" %
                      (tag, method, perr[:40], eerr[:40], gerr[:40]), flush=True)
    perpose.close()

    # ---------------- summary per method --------------------------------
    rows = []
    for method in PROTEIN_METHODS + DNA_METHODS + PYDOCK_ONLY:
        sysname = A[method]
        rec = []
        with open(os.path.join(HERE, "orig_triple_perpose.csv")) as f:
            for r in csv.DictReader(f):
                if r["method"] == method and "ERR" not in r["python"] \
                   and r["python"] != "":
                    rec.append((float(r["python"]), float(r["exact"]),
                                float(r["grid"]), float(r["diff_exact_py"]),
                                float(r["diff_grid_py"])))
        if not rec:
            continue
        py = [x[0] for x in rec]; ex = [x[1] for x in rec]; gr = [x[2] for x in rec]
        dex = [x[3] for x in rec]; dgr = [x[4] for x in rec]
        valid = [i for i in range(len(py)) if abs(py[i]) > 1e-3 and abs(gr[i]) > 1e-3]
        rels = [dgr[i]/abs(py[i]) for i in valid] if valid else []
        rmsd = math.sqrt(sum(d*d for d in dgr)/len(dgr))
        sp_gr = spearman(py, gr)
        sp_ex = spearman(py, ex) if len(set(py)) > 1 else float("nan")
        rows.append({
            "method": method, "system": sysname, "n": len(rec),
            "exact_max_abs_diff": max(dex), "exact_median_abs_diff": statistics.median(dex),
            "grid_max_abs_diff": max(dgr), "grid_median_abs_diff": statistics.median(dgr),
            "grid_max_rel_dev": max(rels) if rels else float("nan"),
            "grid_rmsd": rmsd, "spearman_py_grid": sp_gr, "spearman_py_exact": sp_ex,
        })
        print("\n%s (%s, n=%d): exact-py max %.2e  grid-py max %.3f median %.3f  "
              "grid max_rel %.4f  sp(py,grid)=%.5f" %
              (method, sysname, len(rec), max(dex), max(dgr),
               statistics.median(dgr), max(rels) if rels else float("nan"), sp_gr))
    with open(os.path.join(HERE, "orig_triple_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("\n===== done in %.1fs: exact-vs-original <1e-5 on %d/%d items =====" %
          (time.time() - t_start, n_ok, n_items))
    print("outputs: orig_triple_perpose.csv / orig_triple_summary.csv")


if __name__ == "__main__":
    main()
