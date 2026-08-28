import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})
OUT = "/tmp/lktest/figures"
os.makedirs(OUT, exist_ok=True)

# ---------- Fig A: 评分矩阵热图 ----------
rows = list(csv.DictReader(open("/tmp/lktest/scores/matrix.csv")))
cases = ["1azp", "2oob", "p53DNA", "AbLyso", "AbHIVpep", "RBD_ACE2"]
case_labels = {
    "1azp": "1AZP\n(Rec-DNA)", "2oob": "2OOB\n(Rec-Rec)",
    "p53DNA": "1DIZ\np53-DNA", "AbLyso": "1VFB\nAb-Lysozyme",
    "AbHIVpep": "1DQJ\nAb-HIV peptide", "RBD_ACE2": "6M0J\nRBD-hACE2",
}
methods = ["dfire", "dfire2", "dna", "ddna", "pydock", "cpydock", "sd", "vdw", "mj3h", "pisa", "sipper", "tobi"]
data = np.full((len(methods), len(cases)), np.nan)
for r in rows:
    if r["case"] in cases and r["method"] in methods:
        try:
            data[methods.index(r["method"]), cases.index(r["case"])] = float(r["score"])
        except ValueError:
            pass

fig, ax = plt.subplots(figsize=(9.5, 6.2))
masked = np.ma.masked_invalid(data)
vmin, vmax = np.nanpercentile(data, 2), np.nanpercentile(data, 98)
norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=max(vmax, 1.0))
im = ax.imshow(masked, cmap="RdBu_r", norm=norm, aspect="auto")
ax.set_xticks(range(len(cases)), [case_labels[c] for c in cases], fontsize=8.5)
ax.set_yticks(range(len(methods)), methods, fontsize=9.5)
for i in range(len(methods)):
    for j in range(len(cases)):
        v = data[i, j]
        if np.isnan(v):
            ax.text(j, i, "n/a", ha="center", va="center", fontsize=7.5, color="#999")
        else:
            txt = f"{v:.2f}" if abs(v) < 1000 else f"{v:.1e}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.3,
                    color="white" if abs(v) > abs(vmax)*0.45 else "#222")
cb = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
cb.set_label("Score (lower = more favorable, convention-dependent)", fontsize=8.5)
ax.set_title("LKlight Scoring-Function Matrix: 12 methods x 6 biological complexes\n"
             "(native pose, score module, Rust engine v1.0.2+dfire2-fix)", fontsize=10.5, pad=12)
fig.tight_layout()
fig.savefig(f"{OUT}/figA_score_matrix.png", bbox_inches="tight")
plt.close(fig)

# ---------- Fig B: 参数扫描（g / steps vs 耗时） ----------
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
colors = {"pydock": "#2b6cb0", "dfire": "#c05621"}
import csv as _csv
with open("/tmp/lktest/params/scan.csv") as f:
    lines = [l.strip().split(",") for l in f if l.strip()]
# 跳过表头；行: [method, 'g=25', 'steps=100', '1.865']
for m in ["pydock", "dfire"]:
    gxs, gys, sxs, sys_ = [], [], [], []
    for p in lines[1:]:
        meth, gk, sk, tv = p[0], p[1], p[2], float(p[3])
        if meth != m:
            continue
        if sk == "steps=100":
            gxs.append(int(gk.split("=")[1])); gys.append(tv)
        if gk == "g=200":
            sxs.append(int(sk.split("=")[1])); sys_.append(tv)
    gxs, gys = zip(*sorted(zip(gxs, gys)))
    sxs, sys_ = zip(*sorted(zip(sxs, sys_)))
    axes[0].plot(gxs, gys, "o-", color=colors[m], label=m)
    axes[1].plot(sxs, sys_, "o-", color=colors[m], label=m)
axes[0].set_xlabel("glowworms (g), steps=100"); axes[0].set_ylabel("runtime (s)")
axes[0].set_title("Scaling with swarm size (1AZP)", fontsize=10)
axes[1].set_xlabel("steps, g=200"); axes[1].set_title("Scaling with optimization steps", fontsize=10)
for a in axes:
    a.legend(frameon=False, fontsize=9)
    a.grid(alpha=0.25, ls=":")
fig.suptitle("LKlight parameter sweep: runtime scales linearly with glowworms and steps", fontsize=10.5)
fig.tight_layout()
fig.savefig(f"{OUT}/figB_param_scan.png", bbox_inches="tight")
plt.close(fig)

# ---------- Fig C: 场景耗时对比（docking run, g=200, steps=100） ----------
scen = {
    "p53-DNA\n(1DIZ)": {"dna": 6.15, "ddna": 2.23, "cpydock": 9.77, "dfire": 0.07},
    "Ab-Lysozyme\n(1VFB)": {"pydock": 24.87, "cpydock": 30.79, "dfire": 8.84, "sd": 18.72, "vdw": 7.93},
    "Ab-HIVpep\n(1DQJ)": {"pydock": 38.89, "cpydock": 63.00, "dfire": 42.37, "sd": 55.65, "vdw": 32.68},
    "RBD-hACE2\n(6M0J)": {"pydock": 120.20, "cpydock": 134.28, "dfire": 66.95, "sd": 38.89, "vdw": 34.04},
}
all_methods = ["pydock", "cpydock", "dfire", "sd", "vdw", "dna", "ddna"]
mc = {"pydock": "#2b6cb0", "cpydock": "#805ad5", "dfire": "#c05621", "sd": "#38a169", "vdw": "#d69e2e", "dna": "#3182ce", "ddna": "#9c4221"}
fig, ax = plt.subplots(figsize=(10.5, 4.6))
x = np.arange(len(scen)); width = 0.11
for k, m in enumerate(all_methods):
    vals = [scen[s].get(m, np.nan) for s in scen]
    pos = x + (k - len(all_methods)/2 + 0.5) * width
    bars = ax.bar(pos, np.nan_to_num(vals, nan=0), width*0.92, label=m, color=mc[m])
    for b, v in zip(bars, vals):
        if np.isnan(v):
            b.set_alpha(0.12)
        elif v > 0:
            ax.text(b.get_x()+b.get_width()/2, v+1.5, f"{v:.1f}", ha="center", fontsize=6.8)
ax.set_xticks(x, list(scen.keys()), fontsize=9)
ax.set_ylabel("docking runtime (s), 1 swarm, g=200, 100 steps")
ax.set_title("Multi-scenario docking runtime: 4 biomedical complexes, 7 scoring functions", fontsize=10.5)
ax.legend(frameon=False, fontsize=8.5, ncol=7, loc="upper left")
ax.grid(axis="y", alpha=0.25, ls=":")
# 标注未测组合为浅色
ax.text(0.99, 0.95, "pale bar = combination not run (out-of-scope scoring)", transform=ax.transAxes,
        ha="right", fontsize=7.5, color="#777")
fig.tight_layout()
fig.savefig(f"{OUT}/figC_scenarios.png", bbox_inches="tight")
plt.close(fig)

# ---------- Fig D: 优化收敛（pydock 1AZP, 最终步最优分） ----------
# 从 base 模块 run 收集不同 steps 的 best score
import re
best_by_steps = {}
for st in [10, 25, 50, 100, 200]:
    d = f"/tmp/lktest/params/pydock_s{st}/swarm_0"
    # 最后一个 gso 文件
    import glob
    fs = glob.glob(f"{d}/gso_*.out")
    if not fs:
        continue
    def stepno(p):
        return int(re.search(r"gso_(\d+)\.out", p).group(1))
    last = max(fs, key=stepno)
    scores = [float(l.split()[-1]) for l in open(last) if l.strip() and not l.startswith("#")]
    if scores:
        best_by_steps[st] = min(scores)
if best_by_steps:
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    xs = sorted(best_by_steps)
    ys = [best_by_steps[s] for s in xs]
    ax.plot(xs, ys, "o-", color="#2b6cb0")
    for xx, yy in zip(xs, ys):
        ax.annotate(f"{yy:.1f}", (xx, yy), textcoords="offset points", xytext=(0, -14), ha="center", fontsize=8)
    ax.set_xlabel("GSO optimization steps (1AZP pydock, g=200)")
    ax.set_ylabel("best interface score (final step)")
    ax.set_title("Convergence: best pydock score improves with more GSO steps", fontsize=10.5)
    ax.grid(alpha=0.25, ls=":")
    fig.tight_layout()
    fig.savefig(f"{OUT}/figD_convergence.png", bbox_inches="tight")
    plt.close(fig)
    print("convergence:", best_by_steps)

# ---------- Fig E: 收敛轨迹（单次 100 步 run, 每个输出步的最优分） ----------
import glob
d = "/tmp/lktest/modules/base/swarm_0"
fs = glob.glob(f"{d}/gso_*.out")
steps_scores = []
for p in fs:
    step = int(re.search(r"gso_(\d+)\.out", p).group(1))
    scores = [float(l.split()[-1]) for l in open(p) if l.strip() and not l.startswith("#")]
    if scores:
        steps_scores.append((step, min(scores)))
steps_scores.sort()
if len(steps_scores) >= 3:
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    xs = [s for s, _ in steps_scores]
    ys = [v for _, v in steps_scores]
    ax.plot(xs, ys, "o-", color="#805ad5", ms=4)
    ax.set_xlabel("GSO step (1AZP pydock, 1 swarm, g=200)")
    ax.set_ylabel("best swarm score at checkpoint")
    ax.set_title("Within-run convergence trajectory (checkpoint files gso_N.out)", fontsize=10.5)
    ax.grid(alpha=0.25, ls=":")
    fig.tight_layout()
    fig.savefig(f"{OUT}/figE_convergence_traj.png", bbox_inches="tight")
    plt.close(fig)
    print("traj:", steps_scores[:5], "...", steps_scores[-3:])

print("figures saved to", OUT)
