#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 7 (EN + CN variants): Table 4 visualisation — one-batch 1,000x1,000
engine wall-clocks on the RTX 5090 host (v1.2.0) and the Mac mini M4 host
(CPU engines v1.2.0, Metal backend v1.2.1).

Unified engine naming (manuscript convention):
  - original engine  -> LKlight (reference) / LKlight（参考引擎）
  - CPU accelerated  -> LKlight-grid (CPU build) / LKlight-grid（CPU 版本）
  - GPU accelerated  -> GPU build (LKlight-GPU) / GPU 版本 (LKlight-GPU),
                        with the backend annotated separately: CUDA / Metal.

Outputs:
  figures/fig7_table4_metal.png      (EN labels)
  figures/fig7_table4_metal_cn.png   (CN labels)

Sources (archived): functest/ht1000_summary_5090.csv,
                    functest/ht1000_regression_summary.csv (Mac rows),
                    functest/ht1000_mac/ht1000_summary_mac.csv.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(_HERE)
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "Arial",
    "axes.grid": True, "grid.color": "#d9d9d9", "grid.alpha": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
})

LK_BLUE = "#0072B2"; LK_ORANGE = "#E69F00"; LK_RED = "#D55E00"; LGREY = "#d9d9d9"

def cjk_font():
    for cand in ("Heiti SC", "PingFang SC", "Hiragino Sans GB", "Arial Unicode MS"):
        if any(cand.lower() == f.name.lower() for f in fm.fontManager.ttflist):
            return cand
    return None

# (x label key, wallclock_s) — labels resolved per language
ROWS = [
    ("ref", 128.41), ("grid", 37.01), ("gpu_cuda", 4.05),
    ("ref", 219.29), ("grid", 47.10), ("gpu_metal", 6.22),
]

TEXT = {
    "en": {
        "x": {"ref": "LKlight\n(reference)", "grid": "LKlight-grid\n(CPU)",
              "gpu": "GPU build\n(LKlight-GPU)"},
        "x_backend": {"gpu_cuda": "backend: CUDA", "gpu_metal": "backend: Metal"},
        "legend": ["LKlight (reference)", "LKlight-grid (CPU)",
                   "GPU build — CUDA backend",
                   "GPU build — Metal backend"],
        "backend": {"gpu_cuda": "backend: CUDA", "gpu_metal": "backend: Metal"},
        "hosts": [("Windows + RTX 5090 (v1.2.0)", "GPU: 9.1× vs grid · 31.7× vs LKlight"),
                  ("macOS Mac mini M4, Apple GPU (v1.2.1)", "GPU: 7.6× vs grid · 35.3× vs LKlight")],
        "ylabel": "wall-clock, 1,000 glowworms × 1,000 steps (s, log)",
        "out": "figures/fig7_table4_metal.png",
    },
    "cn": {
        "x": {"ref": "LKlight\n（参考引擎）", "grid": "LKlight-grid\n（CPU）",
              "gpu": "GPU 版本\n(LKlight-GPU)"},
        "x_backend": {"gpu_cuda": "后端：CUDA", "gpu_metal": "后端：Metal"},
        "legend": ["LKlight（参考引擎）", "LKlight-grid（CPU）",
                   "GPU 版本 LKlight-GPU · CUDA 后端",
                   "GPU 版本 LKlight-GPU · Metal 后端"],
        "backend": {"gpu_cuda": "后端：CUDA", "gpu_metal": "后端：Metal"},
        "hosts": [("Windows + RTX 5090（v1.2.0）", "GPU：相对 grid 9.1× · 相对 LKlight 31.7×"),
                  ("macOS Mac mini M4，Apple GPU（v1.2.1）", "GPU：相对 grid 7.6× · 相对 LKlight 35.3×")],
        "ylabel": "墙钟：1,000 萤火虫 × 1,000 步（s，对数轴）",
        "out": "figures/fig7_table4_metal_cn.png",
    },
}

def draw(lang):
    T = TEXT[lang]
    if lang == "cn":
        fam = cjk_font()
        if fam:
            plt.rcParams["font.family"] = fam
    fig, ax = plt.subplots(figsize=(9.4, 4.9))
    xs = np.arange(len(ROWS))
    for i, (key, s) in enumerate(ROWS):
        if key == "ref":
            color, hatch = LGREY, ""
        elif key == "grid":
            color, hatch = LK_BLUE, ""
        elif key == "gpu_cuda":
            color, hatch = LK_RED, ""
        else:
            color, hatch = LK_ORANGE, "//"
        ax.bar(i, s, width=0.62, color=color, edgecolor="#222222",
               linewidth=0.6, hatch=hatch, zorder=3)
        ax.text(i, s * 1.55, f"{s:g}", ha="center", va="bottom",
                fontsize=10, fontweight="bold", color="#111111")
        # backend annotation lives on the 3rd line of the x tick label (see below)

    ax.axvspan(-0.5, 2.5, color="#f4f4f4", zorder=0)
    ax.axvspan(2.5, 5.5, color="#fdf6e3", zorder=0)
    _T = ax.transAxes
    (h1, s1), (h2, s2) = T["hosts"]
    ax.text(0.166, -0.34, h1, ha="center", va="top", transform=_T,
            fontsize=9, fontweight="bold", color="#333333")
    ax.text(0.166, -0.43, s1, ha="center", va="top", transform=_T,
            fontsize=8.6, style="italic", color="#444444")
    ax.text(0.833, -0.34, h2, ha="center", va="top", transform=_T,
            fontsize=9, fontweight="bold", color="#8a6d1a")
    ax.text(0.833, -0.43, s2, ha="center", va="top", transform=_T,
            fontsize=8.6, style="italic", color="#8a6d1a")

    ax.set_xticks(xs)
    def xlabel(k):
        if k in ("ref", "grid"):
            return T["x"][k]
        lab = T["x"]["gpu"] + "\n" + T["x_backend"][k]
        return lab
    ax.set_xticklabels([xlabel(k) for k, _ in ROWS], fontsize=9.2)
    ax.set_ylabel(T["ylabel"])
    ax.set_yscale("log")
    ax.set_ylim(0.6, 600)
    ax.grid(axis="y", which="major", color="#cccccc", linewidth=0.5, zorder=0)
    ax.grid(axis="y", which="minor", color="#e8e8e8", linewidth=0.35, zorder=0)
    ax.minorticks_on()
    ax.tick_params(axis="y", which="minor", length=2.5, color="#bbbbbb")
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    from matplotlib.patches import Patch
    handles = [Patch(facecolor=LGREY, edgecolor="#222222", label=T["legend"][0]),
               Patch(facecolor=LK_BLUE, edgecolor="#222222", label=T["legend"][1]),
               Patch(facecolor=LK_RED, edgecolor="#222222", label=T["legend"][2]),
               Patch(facecolor=LK_ORANGE, edgecolor="#222222", hatch="//", label=T["legend"][3])]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=False, fontsize=8.8, ncol=1)

    os.makedirs("figures", exist_ok=True)
    os.makedirs("figures/pdf", exist_ok=True)
    os.makedirs("figures/tif", exist_ok=True)
    fig.subplots_adjust(bottom=0.44, right=0.70)
    base = T["out"][:-4]
    fig.savefig(f"{base}.png", dpi=300)
    fig.savefig(f"figures/pdf/{os.path.basename(base)}.pdf", bbox_inches="tight")
    fig.savefig(f"figures/tif/{os.path.basename(base)}.tif", bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("wrote", os.path.abspath(f"{base}.png"), "+ pdf/tif")

draw("en")
draw("cn")
