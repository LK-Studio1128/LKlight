#!/usr/bin/env python3
"""出版级图表统一生成脚本（v3 — 终版）。"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "axes.titleweight": "bold",
        "axes.labelweight": "regular",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.frameon": True,
        "legend.fontsize": 11,
        "legend.framealpha": 0.95,
        "legend.edgecolor": "#cccccc",
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.2,
    }
)

C_LK = "#1f4e79"
C_PY = "#2e7d32"
C_RUST = "#c62828"

PAPER = Path("/Users/luoxiaowen/Desktop/LKDock/LKlight论文")
DATA = PAPER / "verify_equivalence"
OUT = PAPER / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def load_30case() -> pd.DataFrame:
    return pd.read_csv(DATA / "bm5_30case" / "summary.tsv", sep="\t")


def load_system_tests() -> dict:
    with (DATA / "system_tests" / "report_data.json").open(encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# 图 2：30 例成功率曲线（终版）
# ---------------------------------------------------------------------------
def plot_fig2_30case(df: pd.DataFrame) -> Path:
    cols = ["top1", "top5", "top10", "top20", "top50", "top100"]
    top_ns = [1, 5, 10, 20, 50, 100]
    lk = df[df.engine == "lk"].set_index("case")[cols].astype(float)
    py = df[df.engine == "py"].set_index("case")[cols].astype(float)

    lk_mean = lk.mean().values
    py_mean = py.mean().values
    # 95% CI（标准误 × 1.96）
    lk_ci = 1.96 * lk.std(ddof=1).values / np.sqrt(len(lk))
    py_ci = 1.96 * py.std(ddof=1).values / np.sqrt(len(py))

    fig, ax = plt.subplots(figsize=(7.0, 4.6))

    ax.fill_between(top_ns, np.maximum(lk_mean - lk_ci, 0), lk_mean + lk_ci,
                    color=C_LK, alpha=0.18, zorder=1)
    ax.fill_between(top_ns, np.maximum(py_mean - py_ci, 0), py_mean + py_ci,
                    color=C_PY, alpha=0.18, zorder=1)
    ax.plot(top_ns, lk_mean, "-o", color=C_LK, lw=2.4, ms=10,
            label="LKlight (this work)", zorder=4, mec="white", mew=1.2)
    ax.plot(top_ns, py_mean, "-s", color=C_PY, lw=2.4, ms=10,
            label="Python LightDock 0.9.4", zorder=4, mec="white", mew=1.2)

    diff = np.abs(lk_mean - py_mean)
    max_d = diff.max()
    # 注释放在左下角
    ax.text(0.02, 0.10, f"max |Δ| = {max_d:.3f}",
            transform=ax.transAxes,
            fontsize=12, color="#444",
            bbox=dict(boxstyle="round,pad=0.4", fc="#fffacd",
                      ec="#b8a064", lw=0.7))

    ax.set_xscale("log")
    ax.set_xticks(top_ns)
    ax.set_xticklabels([str(n) for n in top_ns])
    ax.set_xlim(0.85, 120)
    ax.set_ylim(0, 0.085)
    ax.set_xlabel("Top-N models")
    ax.set_ylabel("Mean success rate")
    ax.set_title("Equivalence on 30 BM5 complexes  (fastdfire, 10 swarms × 50 glowworms)",
                 pad=14)
    ax.legend(loc="upper right", bbox_to_anchor=(0.99, 0.97))

    fig.text(
        0.5, -0.03,
        "Shaded bands = 95% CI over 30 cases.   McNemar p = 0.34,   Wilcoxon p = 0.12 (top-100).",
        ha="center", fontsize=10, color="#555",
    )

    p = OUT / "fig2_30case_success.png"
    fig.savefig(p)
    plt.close(fig)
    return p


# ---------------------------------------------------------------------------
# 图 3：系统测试三面板（终版）
# ---------------------------------------------------------------------------
def _short_tier_label(t: dict) -> str:
    return f"{t['swarms']}×{t['glowworms']}\n×{t['steps']}st"


def plot_fig3_system_tests(payload: dict) -> Path:
    tiers = payload["tiers"]
    scoring = payload["scoring"]
    scenes = payload["scenes"]

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4),
                             gridspec_kw={"width_ratios": [1.1, 1.7, 1.4]})
    fig.subplots_adjust(wspace=0.32, left=0.06, right=0.98, top=0.88, bottom=0.13)

    # -------------------- Panel A: 档位 --------------------
    ax = axes[0]
    labels = [_short_tier_label(t) for t in tiers]
    top5 = [t["top5"] for t in tiers]
    elapsed = [t["elapsed"] for t in tiers]
    x = np.arange(len(tiers))

    bar_colors = ["#bbdefb", "#64b5f6", "#1976d2", "#0d47a1", "#001f4d"]
    bars = ax.bar(x, top5, width=0.55, color=bar_colors,
                  edgecolor="white", linewidth=1.2, zorder=2)
    for b, v in zip(bars, top5):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.022, f"{v:.2f}",
                ha="center", va="bottom", fontsize=11, color="#222", zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, 0.95)
    ax.set_ylabel("top-5 success rate")
    ax.set_title("(a) Parameter tiers  (2X9A, fastdfire)")
    ax.grid(axis="y", zorder=1)

    ax2 = ax.twinx()
    ax2.plot(x, elapsed, "D-", color=C_RUST, lw=2.0, ms=8,
             label="wall-clock", zorder=3, mec="white", mew=0.8)
    for xi, e in zip(x, elapsed):
        ax2.text(xi, e + 1.6, f"{e:.1f}s", ha="center", va="bottom",
                 fontsize=9.5, color=C_RUST)
    ax2.set_ylabel("wall-clock (s)", color=C_RUST)
    ax2.tick_params(axis="y", labelcolor=C_RUST)
    ax2.spines["right"].set_visible(True)
    # 显式 ticks 保证不溢出
    ax2.set_ylim(0, 50)
    ax2.set_yticks([0, 10, 20, 30, 40, 50])

    # 注释（放在左上方空白处，不挡 bar 顶端数字与右侧时间线）
    ax.text(0.05, 0.92,
            "All 5 tiers identical:\nparameter expansion\ndoesn't add new\nsurface points",
            transform=ax.transAxes,
            fontsize=10, color="#222", va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#fffacd",
                      ec="#b8a064", lw=0.7), zorder=6)

    # -------------------- Panel B: 评分函数 --------------------
    ax = axes[1]
    order = sorted(scoring, key=lambda s: -s["top5"])
    names = [s["name"] for s in order]
    top5_vals = [s["top5"] for s in order]
    elapsed_vals = [s["elapsed"] for s in order]
    y = np.arange(len(order))

    cmap = mpl.colormaps.get_cmap("viridis")
    colors = [cmap(0.10 + 0.70 * (i / (len(order) - 1))) for i in range(len(order))]
    bars = ax.barh(y, top5_vals, color=colors, edgecolor="white", lw=1.0, zorder=2)
    for b, v in zip(bars, top5_vals):
        ax.text(v + 0.015, b.get_y() + b.get_height() / 2, f"{v:.2f}",
                va="center", fontsize=10.5, color="#222")
    # 在条形末端右侧标注耗时
    for yi, (b, e) in enumerate(zip(bars, elapsed_vals)):
        # 用一个红钻 + 时间文本放在条形右侧
        ax.text(0.82, yi, f"◆ {e:.1f}s", va="center", ha="left",
                fontsize=10, color=C_RUST, fontweight="medium")
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.15)  # 给右侧时间文字留空间
    ax.set_xlabel("top-5 success rate")
    ax.set_title("(b) Scoring functions  (2X9A)")

    # -------------------- Panel C: 场景气泡 --------------------
    ax = axes[2]
    atoms = np.array([s["atoms"] for s in scenes])
    elapseds = np.array([s["elapsed"] for s in scenes])
    top5s = np.array([s["top5"] for s in scenes])
    labels = [s["name"].split("_")[0] for s in scenes]
    sizes_full = [s["name"].split("_", 1)[1] if "_" in s["name"] else "" for s in scenes]

    type_color = {
        "2X9A": "#1976d2",
        "2GAF": "#7b1fa2",
        "3BIW": "#c62828",
        "1AZP": "#00838f",
        "3MXW": "#ef6c00",
    }
    colors_b = [type_color.get(l, "#666") for l in labels]
    sizes = 250 + 1400 * top5s

    # 标签错开方案：把右上方重叠的三个气泡位置微调
    # 实际数据：2GAF 6120at/101s, 3MXW 4569at/46s, 3BIW 5775at/78s
    # 在 log-log 下它们挤在一起 → 沿用真实数据，但标签加偏移并加边框
    label_offsets = {
        "2GAF": (10, 12),
        "3BIW": (10, -16),
        "3MXW": (-90, -22),
    }

    for x, y, sz, c, lbl, sf, t5 in zip(atoms, elapseds, sizes, colors_b, labels, sizes_full, top5s):
        ax.scatter(x, y, s=sz, c=c, alpha=0.55, edgecolor=c, linewidth=2, zorder=3)
        ox, oy = label_offsets.get(lbl, (8, 6))
        ax.annotate(f"{lbl}  (t5={t5:.2f})", (x, y),
                    xytext=(ox, oy), textcoords="offset points",
                    fontsize=10, color="#222",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white",
                              ec=c, lw=0.8, alpha=0.9))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Complex size (atoms)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title("(c) Application scenarios")
    ax.grid(True, which="both", ls="--", alpha=0.3)

    # O(N²) 趋势线
    ref_idx = int(np.argmin(np.abs(atoms - 1500)))
    ref = (atoms[ref_idx], elapseds[ref_idx])
    xs = np.array([1000, 7000])
    ys = ref[1] * (xs ** 2) / (ref[0] ** 2)
    ax.plot(xs, ys, "--", color="#888", lw=1.4, alpha=0.7)
    ax.text(2900, 15, "O(N²) trend", color="#888", fontsize=10, style="italic")

    # 图例：场景类型
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1976d2',
               markersize=10, label='small PPI (2X9A)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#7b1fa2',
               markersize=10, label='mid PPI (2GAF)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#c62828',
               markersize=10, label='large PPI (3BIW)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#00838f',
               markersize=10, label='protein-DNA (1AZP)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#ef6c00',
               markersize=10, label='antibody-antigen (3MXW)'),
    ]
    ax.legend(handles=legend_elems, loc="lower right", fontsize=9, framealpha=0.95,
              bbox_to_anchor=(0.99, 0.02))

    fig.suptitle(
        "Figure 3 — LKlight system tests: parameter tiers, scoring functions, and application scenarios",
        fontsize=14.5, fontweight="bold", y=1.05,
    )

    p = OUT / "fig3_system_tests.png"
    fig.savefig(p)
    plt.close(fig)
    return p


# ---------------------------------------------------------------------------
# 单图版（备用）
# ---------------------------------------------------------------------------
def plot_tiers_single(payload: dict) -> Path:
    tiers = payload["tiers"]
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    labels = [_short_tier_label(t) for t in tiers]
    top5 = [t["top5"] for t in tiers]
    x = np.arange(len(tiers))
    ax.bar(x, top5, color="#1976d2", edgecolor="white", lw=1)
    for xi, v in zip(x, top5):
        ax.text(xi, v + 0.02, f"{v:.2f}", ha="center", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, 0.95)
    ax.set_ylabel("top-5 success rate")
    ax.set_title("Tier sweep — identical under fixed surface points")
    p = OUT / "fig3a_tiers.png"
    fig.savefig(p)
    plt.close(fig)
    return p


def plot_scoring_single(payload: dict) -> Path:
    scoring = payload["scoring"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    order = sorted(scoring, key=lambda s: -s["top5"])
    names = [s["name"] for s in order]
    top5_vals = [s["top5"] for s in order]
    elapsed_vals = [s["elapsed"] for s in order]
    y = np.arange(len(order))
    cmap = mpl.colormaps.get_cmap("viridis")
    colors = [cmap(0.10 + 0.70 * (i / (len(order) - 1))) for i in range(len(order))]
    ax.barh(y, top5_vals, color=colors, edgecolor="white")
    for yi, v in zip(y, top5_vals):
        ax.text(v + 0.012, yi, f"{v:.2f}", va="center", fontsize=10)
    for yi, e in zip(y, elapsed_vals):
        ax.text(0.83, yi, f"◆ {e:.1f}s", va="center", ha="left",
                fontsize=10, color=C_RUST)
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("top-5 success rate")
    ax.set_title("Scoring function sweep (2X9A)")
    p = OUT / "fig3b_scoring.png"
    fig.savefig(p)
    plt.close(fig)
    return p


def plot_scenes_single(payload: dict) -> Path:
    scenes = payload["scenes"]
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    atoms = [s["atoms"] for s in scenes]
    elapseds = [s["elapsed"] for s in scenes]
    top5s = [s["top5"] for s in scenes]
    labels = [s["name"].split("_")[0] for s in scenes]
    sizes = 250 + 1400 * np.array(top5s)
    type_color = {
        "2X9A": "#1976d2", "2GAF": "#7b1fa2", "3BIW": "#c62828",
        "1AZP": "#00838f", "3MXW": "#ef6c00",
    }
    label_offsets = {
        "2GAF": (10, 12),
        "3BIW": (10, -16),
        "3MXW": (-90, -22),
    }
    for x, y, sz, lbl, t5 in zip(atoms, elapseds, sizes, labels, top5s):
        ax.scatter(x, y, s=sz, alpha=0.55,
                   edgecolor=type_color.get(lbl, "#444"), linewidth=2)
        ox, oy = label_offsets.get(lbl, (8, 6))
        ax.annotate(f"{lbl} (t5={t5:.2f})", (x, y),
                    xytext=(ox, oy), textcoords="offset points",
                    fontsize=10,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white",
                              ec=type_color.get(lbl, "#444"), lw=0.8, alpha=0.9))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Complex size (atoms)")
    ax.set_ylabel("Wall-clock (s)")
    ax.set_title("Application scenarios — size × time × success")
    ax.grid(True, which="both", ls="--", alpha=0.3)
    ref_idx = int(np.argmin(np.abs(np.array(atoms) - 1500)))
    ref = (atoms[ref_idx], elapseds[ref_idx])
    xs = np.array([1000, 7000])
    ys = ref[1] * (xs ** 2) / (ref[0] ** 2)
    ax.plot(xs, ys, "--", color="#888", lw=1.4, alpha=0.7)
    ax.text(2900, 15, "O(N²) trend", color="#888", fontsize=10, style="italic")
    p = OUT / "fig3c_scenes.png"
    fig.savefig(p)
    plt.close(fig)
    return p


def main() -> None:
    df30 = load_30case()
    payload = load_system_tests()

    paths = [
        plot_fig2_30case(df30),
        plot_fig3_system_tests(payload),
        plot_tiers_single(payload),
        plot_scoring_single(payload),
        plot_scenes_single(payload),
    ]
    for p in paths:
        print(f"wrote {p}  ({p.stat().st_size / 1024:.1f} KiB)")


if __name__ == "__main__":
    main()