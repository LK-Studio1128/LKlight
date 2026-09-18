#!/usr/bin/env python3
"""Server-portable equivalence harness: base(exact) vs CPU-grid.

Identical logic/metrics/seed to equiv_harness.py (Mac), but BASE / GRID / REC / LIG
are taken from environment variables so the same script runs on the Linux server
without editing ROOT paths.  Run on the server:

    BASE=<exact bin> GRID=<grid bin> \
    REC=<rec.pdb>    LIG=<lig.pdb>   \
    python3 equiv_harness_server.py --nposes 30 --out equiv_linux

Outputs <out>_summary.csv and <out>_perpose.csv in the script's directory.
"""
import argparse, csv, math, random, statistics, subprocess, time, os, sys

BASE   = os.environ.get("BASE")     # all-pairs exact engine
GRID   = os.environ.get("GRID")     # CPU-grid engine (or Option-A build)
REC    = os.environ.get("REC")
LIG    = os.environ.get("LIG")

METHODS = ["dfire","dfire2","dna","mj3h","pydock","cpydock","sd","vdw",
           "pisa","sipper","tobi","ddna"]

def normq(qw,qx,qy,qz):
    n = math.sqrt(qw*qw+qx*qx+qy*qy+qz*qz) or 1.0
    return qw/n, qx/n, qy/n, qz/n

def quat_from_axis(axis, deg):
    r = math.radians(deg); s = math.sin(r/2); c = math.cos(r/2)
    ax,ay,az = axis
    return normq(c, ax*s, ay*s, az*s)

def make_decoys(n):
    """native (0,0,0, identity) + n rigid-body decoys: random translation
    (0.5..12 A, biased so many stay within the ~10 A window plus some that
    probe the far-field 10-30 A band) and rotation up to 60 deg."""
    decoys = [(0.0,0.0,0.0,1.0,0.0,0.0,0.0)]
    rng = random.Random(20260904)
    for _ in range(n):
        r = rng.uniform(0,1)
        mag = rng.uniform(2.0, 8.0) if r < 0.6 else rng.uniform(10.0, 26.0)
        th = rng.uniform(0, 2*math.pi); ph = math.acos(rng.uniform(-1,1))
        tx,ty,tz = mag*math.sin(ph)*math.cos(th), mag*math.sin(ph)*math.sin(th), mag*math.cos(ph)
        deg = rng.uniform(0, 60)
        ax = [rng.uniform(-1,1), rng.uniform(-1,1), rng.uniform(-1,1)]
        L = math.sqrt(sum(a*a for a in ax)) or 1.0
        ax = [a/L for a in ax]
        qw,qx,qy,qz = quat_from_axis(ax, deg)
        decoys.append((tx,ty,tz,qw,qx,qy,qz))
    return decoys

def run(binpath, method, tx,ty,tz,qw,qx,qy,qz):
    args = [binpath,"score",REC,LIG,method,
            "--tx",f"{tx:.6f}","--ty",f"{ty:.6f}","--tz",f"{tz:.6f}",
            "--qw",f"{qw:.6f}","--qx",f"{qx:.6f}","--qy",f"{qy:.6f}","--qz",f"{qz:.6f}"]
    t0 = time.perf_counter()
    p = subprocess.run(args, capture_output=True, text=True)
    dt = time.perf_counter()-t0
    out = (p.stdout or "")+(p.stderr or "")
    for ln in out.splitlines():
        if ln.startswith("Score ("):
            val = float(ln.split(":")[1].strip())
            return val, dt
    return float("nan"), dt

def pearson(x,y):
    n=len(x); mx=sum(x)/n; my=sum(y)/n
    cov=sum((x[i]-mx)*(y[i]-my) for i in range(n))
    vx=math.sqrt(sum((v-mx)**2 for v in x)); vy=math.sqrt(sum((v-my)**2 for v in y))
    return cov/(vx*vy) if vx*vy else float("nan")

def spearman(x,y):
    def rk(s):
        idx=sorted(range(len(s)), key=lambda i:s[i])
        r=[0]*len(s); rank=0; i=0
        while i<len(s):
            j=i
            while j+1<len(s) and s[idx[j+1]]==s[idx[i]]: j+=1
            a=(i+j)/2
            for k in range(i,j+1): r[idx[k]]=a
            i=j+1
        return r
    rx,ry=rk(x),rk(y); return pearson(rx,ry)

def kendall(x,y):
    c=0;d=0
    for i in range(len(x)):
        for j in range(i+1,len(x)):
            a=(x[i]-x[j])*(y[i]-y[j])
            if a>0:c+=1
            elif a<0:d+=1
    return (c-d)/(c+d) if (c+d) else float("nan")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--nposes",type=int,default=28)
    ap.add_argument("--out",default="equiv_server")
    a=ap.parse_args()
    for name,val in [("BASE",BASE),("GRID",GRID),("REC",REC),("LIG",LIG)]:
        if not val or not os.path.exists(val):
            sys.exit(f"ERROR: env {name} missing or not found: {val!r}")
    decoys=make_decoys(a.nposes)
    here=os.path.dirname(os.path.abspath(__file__))
    rows=[]; perpose=[]
    for m in METHODS:
        E=[];G=[];poses=[]
        for (tx,ty,tz,qw,qx,qy,qz) in decoys:
            e,_=run(BASE,m,tx,ty,tz,qw,qx,qy,qz)
            g,_=run(GRID,m,tx,ty,tz,qw,qx,qy,qz)
            E.append(e);G.append(g)
            poses.append((tx,ty,tz,qw,qx,qy,qz))
            perpose.append({"method":m,"pose":len(E)-1,
                            "tx":round(tx,3),"ty":round(ty,3),"tz":round(tz,3),
                            "exact":e,"grid":g,"diff":g-e})
        tol=1e-6
        valid=[i for i in range(len(E)) if abs(E[i])>tol and abs(G[i])>tol]
        if not valid:
            valid=list(range(len(E)))
        x=[E[i] for i in valid]; y=[G[i] for i in valid]
        devs=[abs(y[k]-x[k])/abs(x[k]) if abs(x[k])>1e-3 else abs(y[k]-x[k]) for k in range(len(x))]
        rmsd=math.sqrt(sum((y[k]-x[k])**2 for k in range(len(x)))/len(x))
        bitid=all(abs(y[k]-x[k])<1e-9 for k in range(len(x)))
        sp=spearman(x,y); pe=pearson(x,y); kt=kendall(x,y)
        maxrel=max(devs) if devs else float("nan")
        if bitid: verdict="bit-identical"
        elif maxrel<0.01 and sp>0.999: verdict="grid-eq (<1%)"
        elif maxrel<0.05 and sp>0.99:  verdict="grid-eq (<5%)"
        else:                          verdict="CHECK"
        rows.append({"method":m,"bit_identical":bitid,"max_rel_dev":maxrel,
                     "rmsd":rmsd,"spearman":sp,"pearson":pe,"kendall":kt,
                     "verdict":verdict,"n_valid":len(x)})
        print(f"{m:8s} {verdict:16s} max|rel|={maxrel:9.4f}  sp={sp:.4f} rmsd={rmsd:.3e}",
              flush=True)
    with open(os.path.join(here,a.out+"_summary.csv"),"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows: w.writerow(r)
    with open(os.path.join(here,a.out+"_perpose.csv"),"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(perpose[0].keys())); w.writeheader()
        for r in perpose: w.writerow(r)
    print(f"\nWrote {a.out}_summary.csv and {a.out}_perpose.csv in {here}")

if __name__=="__main__":
    main()
