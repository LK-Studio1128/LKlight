#!/usr/bin/env python3
"""GSO steady-state timing + convergence: base(exact) vs CPU-grid, per method.

For each of the 12 scoring functions, run a full GSO optimization (20 glow x
100 steps, fixed seed) through the base all-pairs binary and the grid binary
in a fresh working copy, time it, and parse the best (minimum) trajectory
"Scoring" energy across all glowworm output files as the docking-convergence
score. Reports per method: wall time, best energy, relative delta of best.

The whole-run wall time amortises the one-time field/cell construction, so it
is a true steady-state measure (per-process `score` timing is NOT used).

Usage: python3 gso_harness.py [--steps N] [--out prefix]
"""
import argparse, csv, glob, os, re, shutil, statistics, subprocess, time

HERE   = os.path.dirname(os.path.abspath(__file__))
ROOT   = "/Users/luoxiaowen/Desktop/LKDock/byi"
BASE   = f"{ROOT}/LKlight/target/release/LKlight"
GRID   = f"{ROOT}/LKlight-grid/target/release/LKlight"
REC    = os.path.join(HERE,"1azp_receptor.pdb")
LIG    = os.path.join(HERE,"1azp_ligand.pdb")
BASE_SETUP = os.path.join(HERE,"setup.json")   # 1 swarm x 20 glow, already built
METHODS = ["dfire","dfire2","dna","mj3h","pydock","cpydock","sd","vdw",
           "pisa","sipper","tobi","ddna"]

def best_score(workdir):
    """min trajectory Scoring across all glowworm gso_*.out in swarm_0/."""
    pat=re.compile(r"^\s*[^#].*?\s(-?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*$")
    best=None
    files=glob.glob(os.path.join(workdir,"swarm_0","gso_*.out"))
    for f in files:
        for line in open(f,errors="ignore"):
            if line.startswith("#") or not line.strip(): continue
            try:
                val=float(line.split()[-1])   # last column = Scoring
            except ValueError:
                continue
            if best is None or val<best: best=val
    return best

def run_one(binpath, method, workdir, steps):
    shutil.rmtree(workdir,ignore_errors=True); os.makedirs(workdir,exist_ok=True)
    shutil.copy(BASE_SETUP, os.path.join(workdir,"setup.json"))
    dat=os.path.join(HERE,"initial_positions_0.dat")
    shutil.copy(dat, os.path.join(workdir,"initial_positions_0.dat"))
    # the run command reads the lightdock_-prefixed receptor/ligand (as produced by `setup`)
    for p in ["lightdock_1azp_receptor.pdb","lightdock_1azp_ligand.pdb"]:
        shutil.copy(os.path.join(HERE,p), os.path.join(workdir,p))
    t0=time.perf_counter()
    p=subprocess.run([binpath,"run",os.path.join(workdir,"setup.json"),
                      os.path.join(workdir,"initial_positions_0.dat"),
                      str(steps),method],cwd=workdir,
                     capture_output=True,text=True)
    dt=time.perf_counter()-t0
    if p.returncode!=0:
        return dt, None, (p.stderr or "")[-200:]
    return dt, best_score(workdir), None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--steps",type=int,default=100)
    ap.add_argument("--out",default="gso_results")
    ap.add_argument("--root",default=HERE,
        help="working-dir root (use /tmp to avoid protected-path bulk deletes)")
    a=ap.parse_args()
    os.makedirs(a.root,exist_ok=True)
    rows=[]
    for m in METHODS:
        wt=os.path.join(a.root,f"_w_{m}")
        gt=os.path.join(a.root,f"_w_{m}_grid")
        tb,b,eb=run_one(BASE,m,wt,a.steps)
        tg,g,eg=run_one(GRID,m,gt,a.steps)
        rel = None
        if b is not None and g is not None and abs(b)>1e-3:
            rel=abs(g-b)/abs(b)
        elif b is not None and g is not None:
            rel=abs(g-b)
        sp=(tb/tg) if tg>0 else None
        rows.append({"method":m,"steps":a.steps,
                     "base_s":round(tb,4),"grid_s":round(tg,4),
                     "speedup": (round(sp,3) if sp else ""),
                     "base_best": (round(b,4) if b is not None else None),
                     "grid_best": (round(g,4) if g is not None else None),
                     "best_rel_delta": (round(rel,6) if rel is not None else None),
                     "err": (eb or eg or "")})
        print(f"{m:8s} base={tb:6.3f}s grid={tg:6.3f}s speedup={sp if sp else float('nan'):5.2f}x "
              f"base_best={b if b is not None else float('nan'):.1f} grid_best={g if g is not None else float('nan'):.1f} "
              f"best_delta={rel if rel is not None else float('nan'):.2e}")
    with open(os.path.join(HERE,a.out+"_summary.csv"),"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows: w.writerow(r)
    print("\nWrote",a.out+"_summary.csv")

if __name__=="__main__":
    main()
