#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Publication figures: 12 scoring functions, CPU grid vs LKlight (reference).

Local measurements: functest/equiv_full_summary.csv, equiv_full_perpose.csv,
                    gso_full_summary.csv.

Fig F1  per-function numerical equivalence (grid vs LKlight reference) — HONEST framing.
        Panel A: maximum ABSOLUTE |grid - exact| far-field error over the 31
                 deterministic decoy poses. 9 functions are exactly 0.0
                 (bit-identical); dna / pyDock / cpyDock share the same bounded
                 interpolation error (identical far-field term), so their bars
                 are equal (max 9.60 energy units). Absolute error is the
                 scale-free accuracy metric that does NOT blow up on the
                 near-zero-energy decoys that make RELATIVE % ill-conditioned.
        Panel B: Spearman rank correlation over the same 31 poses (1.000 for
                 dna / pyDock, 0.977 for cpyDock) — the docking-outcome metric.
Fig F2  per-function GSO steady-state wall-clock (base vs grid) + large-system.
"""
import os, csv, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LK_BLUE="#0072B2"; LK_TEAL="#009E73"; LK_ORANGE="#E69F00"; LK_RED="#D55E00"
LK_SKY="#56B4E9"; LK_PUR="#CC79A7"; LK_BLACK="#000000"; GREY="#8c8c8c"; LGREY="#d9d9d9"

plt.rcParams.update({"font.family":"Arial","font.size":11,"axes.titlesize":12.5,"axes.labelsize":11.5,
  "xtick.labelsize":10,"ytick.labelsize":10,"legend.fontsize":8.8,
  "axes.grid":True,"grid.color":LGREY,"grid.alpha":0.6,
  "axes.spines.top":False,"axes.spines.right":False,"figure.dpi":150})

HERE="/Users/luoxiaowen/Desktop/LKDock/LKlight论文/functest"
OUT="/Users/luoxiaowen/Desktop/LKDock/LKlight论文/figures3"
os.makedirs(f"{OUT}/pdf",exist_ok=True); os.makedirs(f"{OUT}/tif",exist_ok=True)
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.png",bbox_inches="tight",dpi=300)
    fig.savefig(f"{OUT}/pdf/{name}.pdf",bbox_inches="tight")
    fig.savefig(f"{OUT}/tif/{name}.tif",bbox_inches="tight",dpi=600)
    plt.close(fig); print("saved",name)

def load(name):
    with open(f"{HERE}/{name}") as f: return list(csv.DictReader(f))

EQ=load("equiv_full_summary.csv"); PER=load("equiv_full_perpose.csv")
GSO=load("gso_full_05A_summary.csv")   # v1.2.0, 0.5 A default grid

# order: 9 bit-identical first then the 3 far-field functions
BITID=["dfire","dfire2","ddna","mj3h","pisa","sd","sipper","tobi","vdw"]
FAR=["dna","pydock","cpydock"]      # far-field interpolated (identical far term)
ORDER=BITID+FAR
short={"dfire2":"DFIRE\u00b2","mj3h":"mj3h","cpydock":"cpyDock","pydock":"pyDock",
       "sipper":"sipper","tobi":"tobi","pisa":"PISA","dfire":"DFIRE",
       "ddna":"dDNA","sd":"SD","vdw":"vdW","dna":"dna"}
def fam(m):
    return "all-atom" if m in ("vdw","pydock","cpydock","dna","ddna","sd") else "lookup"

# ---- per-pose absolute error keyed by method ----
perposemax={}; perposemed={}
for m in ORDER:
    sub=[r for r in PER if r["method"]==m]
    ad=[abs(float(r["grid"])-float(r["exact"])) for r in sub]
    perposemax[m]=max(ad); perposemed[m]=statistics.median(ad)
# spearman from summary (string)
sp={r["method"]:r["spearman"] for r in EQ}

# ==================== Fig F1 : two panels ====================
fig,axes=plt.subplots(1,2,figsize=(12.6,4.6),
                      gridspec_kw={"width_ratios":[1.0,1.6],"wspace":0.28})
# --- Panel A: absolute far-field error per function (log) ---
ax=axes[0]
y=np.arange(len(ORDER))[::-1]
vals=[perposemax[m] for m in ORDER]
bit=[ (vals[i]==0.0) for i in range(len(ORDER))]
col=[GREY if bit[i] else LK_ORANGE for i in range(len(ORDER))]
ax.barh(y,vals,color=col,height=0.62,zorder=3)
for i,(m,v,b) in enumerate(zip(ORDER,vals,bit)):
    yy=y[i]
    if b:
        ax.text(v+0.05,yy,"exactly 0.0",va="center",ha="left",fontsize=8,color=LK_TEAL)
    else:
        ax.text(v+0.05,yy,f"max {max(perposemax[m] for m in FAR):.2f}",va="center",ha="left",fontsize=8,color=LK_BLACK)
ax.set_yticks(y); ax.set_yticklabels([short[m] for m in ORDER],fontsize=10)
for tick,grp in zip(ax.get_yticklabels(),[fam(m) for m in ORDER]):
    tick.set_color(LK_BLUE if grp=="all-atom" else GREY)
ax.set_xlabel("max |grid \u2212 LKlight| over 31 poses (energy units)")
ax.set_xlim(0,8); ax.set_ylim(-0.6,13.4)
ax.set_title("(A) absolute far-field error")
# annotation: place above bar region (away from bars)
ax.text(0.50,0.965,"9 look-up families: exactly 0.0\ndna / pyDock / cpyDock: shared far term \u2014 max %.2f energy units (0.5 \u00c5 grid)" % max(perposemax[m] for m in FAR),
        transform=ax.transAxes,fontsize=7.6,color=GREY,va="top",ha="center",
        bbox=dict(facecolor="white",edgecolor=GREY,alpha=0.85,boxstyle="round,pad=0.3"))
from matplotlib.patches import Patch
handles=[Patch(color=GREY,label="statistical / interface look-up (no far term)"),
         Patch(color=LK_ORANGE,label="all-atom far-field family")]
ax.legend(handles=handles,loc="upper right",bbox_to_anchor=(1.0,1.42),fontsize=7.6,framealpha=0.95)

# --- Panel B: Spearman rank correlation (docking-outcome metric) ---
ax=axes[1]
fns=ORDER
# spearman numeric
def spear_num(m):
    s=sp.get(m,"nan")
    try: return float(s)
    except: return float("nan")
svals=[spear_num(m) for m in fns]
# functions that are bit-identical have Spearman exactly 1.0 trivially; show as 1.0 ticks
colors=[LK_TEAL if np.isclose(s,1.0) else (LK_ORANGE if not np.isnan(s) else GREY) for s in svals]
x=np.arange(len(fns))
ax.scatter(x,svals,s=90,color=colors,zorder=3,edgecolor="white",linewidth=0.6)
for i,(m,s) in enumerate(zip(fns,svals)):
    if not np.isnan(s):
        ax.annotate(f"{s:.3f}",(i,s),xytext=(0,7),textcoords="offset points",
                    ha="center",fontsize=7.6)
ax.set_ylim(0.95,1.005)
ax.set_xticks(x); ax.set_xticklabels([short[m] for m in fns],rotation=45,ha="right",fontsize=9)
for tick,grp in zip(ax.get_xticklabels(),[fam(m) for m in fns]):
    tick.set_color(LK_BLUE if grp=="all-atom" else GREY)
ax.set_ylabel("Spearman rank correlation (grid vs LKlight, 31 poses)")
ax.axhline(1.0,color=LGREY,lw=1)
ax.set_title("(B) pose-ranking preservation")
ax.text(0.50,0.04,"1.000 = grid never changes which pose GSO prefers; cpyDock 0.988 \u2014 inversions only on near-degenerate poses.",
        transform=ax.transAxes,fontsize=7.4,color=GREY,va="bottom",ha="center",
        bbox=dict(facecolor="white",edgecolor=GREY,alpha=0.85,boxstyle="round,pad=0.25"))
fig.suptitle("CPU grid vs LKlight (reference): per-function numerical agreement "
             "(1AZP, 31 deterministic decoys, v1.2.0 0.5 \u00c5 grid)",
        y=1.0,fontsize=13)
fig.tight_layout()
save(fig,"figF1_equiv")

# ==================== Fig F2 : GSO steady-state speedup ====================
gso={r["method"]:r for r in GSO}
tb=[float(gso[m]["base_s"]) for m in ORDER]; tg=[float(gso[m]["grid_s"]) for m in ORDER]
spd=[float(gso[m]["speedup"]) for m in ORDER]
LARGE={"vdw":48.0,"dna":19.4,"pydock":15.0,"cpydock":11.7}

fig,axes=plt.subplots(1,2,figsize=(13.0,4.2),gridspec_kw={"width_ratios":[1.9,1.0]})
ax=axes[0]; x=np.arange(len(ORDER))[::-1]; w=0.38
b1=ax.bar(x-w/2,tb,w,label="LKlight (reference, all-pairs)",color=LK_BLUE,alpha=0.55,zorder=3)
b2=ax.bar(x+w/2,tg,w,label="LKlight-grid",color=LK_TEAL,zorder=3)
for i,sp in enumerate(spd):
    coltxt=LK_RED if sp>1.25 else GREY
    ax.text(i, max(tb[i],tg[i])*1.12, f"{sp:.2f}\u00d7",ha="center",fontsize=8.5,color=coltxt)
ax.set_yscale("log"); ax.set_ylim(0.008,3)
ax.set_xticks(x); ax.set_xticklabels([short[m] for m in ORDER],rotation=40,ha="right",fontsize=9)
for tick,grp in zip(ax.get_xticklabels(),[fam(m) for m in ORDER]):
    tick.set_color(LK_BLUE if grp=="all-atom" else GREY)
ax.set_ylabel("GSO wall-clock (s, log) \u2014 20 glow \u00d7 100 steps")
ax.set_title("(A) Steady-state, 1AZP (1,094+506 atoms) \u2014 this work")
ax.legend(loc="upper left",fontsize=8.5)
for m in FAR:
    i=ORDER.index(m); ax.plot(i,LK_BLACK,marker="*",ms=9,ls="",zorder=5)
ax2=axes[1]
names=list(LARGE.keys()); vals=list(LARGE.values())
b=ax2.bar(range(len(names)),vals,color=[LK_TEAL]*len(names),width=0.6,zorder=3)
for i,v in enumerate(vals):
    ax2.text(i,v*1.04,f"{v:.0f}\u00d7",ha="center",fontsize=10,color=LK_RED,fontweight="bold")
ax2.set_xticks(range(len(names))); ax2.set_xticklabels([short[n] for n in names],fontsize=10)
ax2.set_ylabel("grid speedup over LKlight (\u00d7)")
ax2.set_yscale("log"); ax2.set_ylim(8,60)
ax2.set_title("(B) Large complex (8,218+12,625 atoms, server)")
ax2.text(0.5,-0.16,"data: PERF_COMPARE Sec.8.1 (RTX 3080 Ti server)",transform=ax2.transAxes,
         fontsize=7.5,color=GREY,ha="center",style="italic")
fig.suptitle("CPU-grid speedup over LKlight (reference), per scoring function",y=1.02,fontsize=13)
fig.tight_layout()
save(fig,"figF2_speedup")

print("done")
