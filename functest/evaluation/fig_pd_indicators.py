#!/usr/bin/env python3
"""Figure for the protein-DNA indicator case study (1AZP, Sac7d-DNA).

Panel A: distribution of ligand RMSD over the pooled decoy sets, ANM on vs
         off, with the CAPRI L-RMSD thresholds marked.
Panel B: best DockQ as a function of the top-N cut-off, for the scoring-based
         ranking (what a run reports) and for the oracle ranking (decoys
         ordered by true L-RMSD, i.e. the sampling ceiling).

Reads the eval_metrics.csv written by eval_run.py.
"""
import sys
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # functest/
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dockq_eval as _de  # noqa: E402
OUT = os.path.join(os.path.dirname(ROOT), "figures3")               # LKlight论文/figures3
plt.rcParams.update({"font.family": "Arial", "font.size": 9,
                     "axes.edgecolor": "#333333", "axes.linewidth": 0.9})


def load(tag):
    rows = []
    with open(os.path.join(ROOT, "e4_dna", tag, "eval_metrics.csv")) as fh:
        for r in csv.DictReader(fh):
            rows.append(dict(lrmsd=float(r["lrmsd"]), fnat=float(r["fnat"]),
                             dockq=float(r["dockq"]), capri=r["capri"],
                             score=float(r["score"])))
    return rows


def curve(rows, ns, key=lambda r: r["score"], reverse=False):
    seq = sorted(rows, key=key, reverse=reverse)
    out = []
    for n in ns:
        sub = seq[:n]
        if not sub:
            continue
        out.append((n, max(s["dockq"] for s in sub),
                    sum(1 for s in sub if s["capri"] != "Incorrect")))
    return out


on, off = load("anm_on"), load("anm_off")
BLUE, ORANGE, GREY = "#2166ac", "#d95f02", "#777777"
fig = plt.figure(figsize=(13.4, 3.9), dpi=300)
axA = fig.add_subplot(1, 3, 1)
axB = fig.add_subplot(1, 3, 2)
axC = fig.add_subplot(1, 3, 3, projection='3d')

# ── Panel A: L-RMSD distributions ────────────────────────────────────────────
bins = np.arange(0, 45, 1.0)
axA.hist([r["lrmsd"] for r in off], bins=bins, color=ORANGE, alpha=0.70,
         label=f"ANM off (n={len(off)})")
axA.hist([r["lrmsd"] for r in on], bins=bins, color=BLUE, alpha=0.62,
         label=f"ANM on (n={len(on)})")
for thr, lab in ((1.0, "1"), (2.0, "2"), (4.0, "4")):
    axA.axvline(thr, color=GREY, ls=":", lw=1.0)
axA.text(4.2, axA.get_ylim()[1] * 0.93, "CAPRI L-RMSD\nthresholds 1/2/4 Å",
         fontsize=7.0, color=GREY)
c_on, c_off = min(on, key=lambda r: r["lrmsd"]), min(off, key=lambda r: r["lrmsd"])
axA.annotate(f"sampling ceiling\nANM on: {c_on['lrmsd']:.1f} Å\nANM off: {c_off['lrmsd']:.1f} Å",
             xy=(c_on["lrmsd"], 0), xytext=(15.5, axA.get_ylim()[1] * 0.44),
             fontsize=7.4, color="#222222",
             arrowprops=dict(arrowstyle="->", color="#222222", lw=0.9))
axA.set_xlabel("L-RMSD of decoy vs bound complex (Å)")
axA.set_ylabel("number of decoys")
axA.set_title("A  Protein–DNA decoy quality (1AZP, Sac7d–DNA;\n"
              "10 swarms × 200 glowworms × 200 steps)", fontsize=9, loc="left")
axA.legend(frameon=False, fontsize=7.8, loc="upper right")

# ── Panel B: best DockQ vs top-N ─────────────────────────────────────────────
ns = [1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
for rows, color, name in ((on, BLUE, "ANM on"), (off, ORANGE, "ANM off")):
    scored = curve(rows, ns)
    oracle = curve(rows, ns, key=lambda r: r["lrmsd"])
    axB.plot([s[0] for s in scored], [s[1] for s in scored], "o-", color=color,
             ms=3.4, lw=1.5, label=f"{name} — scoring-ranked")
    axB.plot([s[0] for s in oracle], [s[1] for s in oracle], "s--", color=color,
             ms=3.0, lw=1.2, alpha=0.65, label=f"{name} — oracle (ceiling)")
axB.axhline(0.23, color=GREY, ls=":", lw=1.0)
axB.text(2000, 0.243, "DockQ of a CAPRI-acceptable model (≈0.23)", fontsize=6.8,
         color=GREY, ha="right")
axB.set_xscale("log")
axB.set_xlabel("top-N cut-off (pooled decoys, log scale)")
axB.set_ylabel("best DockQ within top-N")
axB.set_ylim(0, 0.33)
axB.set_title("B  Ranking vs sampling ceiling", fontsize=9, loc="left")
axB.legend(frameon=False, fontsize=6.6, loc="upper left", borderaxespad=0.2)

# ── Panel C: 3D overlay of the best decoy per configuration ─────────────────
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import csv as _csv

def best_decoy(tag):
    rows = list(_csv.DictReader(open(os.path.join(ROOT, "e4_dna", tag, "eval_metrics.csv"))))
    best = min(rows, key=lambda r: float(r["lrmsd"]))
    return os.path.join(ROOT, "e4_dna", tag, best["decoy"]), float(best["lrmsd"])

import dockq_eval as _de
nat_rec = _de.parse_pdb(os.path.join(ROOT, "e4_dna", "native", "1azp_protein.pdb"))
nat_lig = _de.parse_pdb(os.path.join(ROOT, "e4_dna", "native", "1azp_dna.pdb"))
for atoms, color, size, label in ((nat_rec, "#bbbbbb", 2.5, "native protein"),
                                  (nat_lig, "#222222", 4.0, "native DNA")):
    xs = [a[4][0] for a in atoms]; ys = [a[4][1] for a in atoms]; zs = [a[4][2] for a in atoms]
    axC.scatter(xs, ys, zs, s=size, c=color, alpha=0.55, depthshade=False,
                label=label, linewidths=0)
for tag, color, name in (("anm_off", ORANGE, "best decoy, ANM off"),
                         ("anm_on", BLUE, "best decoy, ANM on")):
    path, lr = best_decoy(tag)
    atoms = _de.parse_pdb(path)
    lig = [a for a in atoms if a[0] != "A"]
    xs = [a[4][0] for a in lig]; ys = [a[4][1] for a in lig]; zs = [a[4][2] for a in lig]
    axC.scatter(xs, ys, zs, s=5.0, c=color, depthshade=False, linewidths=0,
                label=f"{name} (L-RMSD {lr:.1f} Å)")
axC.set_xlabel("x (Å)", fontsize=6.5, labelpad=-4); axC.set_ylabel("y (Å)", fontsize=6.5, labelpad=-4)
axC.set_zlabel("z (Å)", fontsize=6.5, labelpad=-4)
axC.tick_params(labelsize=6, pad=-2)
axC.set_title("C  Best decoy vs bound complex", fontsize=9, loc="left")
axC.legend(frameon=False, fontsize=6.2, loc="upper left", bbox_to_anchor=(0.72, 0.95))
axC.view_init(elev=18, azim=-60)
axC.set_box_aspect((1, 1, 0.8))

for ax in (axA, axB):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

fig.tight_layout()
os.makedirs(OUT, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"fig_pd_indicators.{ext}"), dpi=300,
                bbox_inches="tight", facecolor="white")
print("wrote", os.path.join(OUT, "fig_pd_indicators.png/pdf"))
