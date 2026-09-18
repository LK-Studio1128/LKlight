# LKlight 测试数据归档（FINAL_DATA）

统一整理版：论文《LKlight: a complete, quantitatively benchmarked Rust engine for the
LightDock protein docking protocol》全部测试数据的可复现副本。与论文正文、图 1–3、
表 1–5 一一对应。

**生成日期：** 2026-08-28
**数据版本：** v1.1（对应论文投稿版）
**原始位置：** `LKlight论文/` 各子目录（本目录为快照副本，保证投稿后可独立引用）

---

## 目录总览

| 子目录 | 内容 | 对应论文 |
|---|---|---|
| `01_numeric_validation/` | 48 项数值对比（12 评分函数 × 4 刚体 pose） | §3.1 表 2 |
| `02_bm5_30case/` | 30 例 BM5 CAPRI 成功率（fastdfire） | §3.3 表 4、图 2 |
| `03_system_tests/` | 档位/评分函数/场景三维系统测试 | §3.4 表 5、图 3 |
| `04_win_crossplat/` | Windows x86-64 跨平台字节级一致性 | §3.1 表 2、§1.4 |
| `05_performance/` | 性能基准原始计时 | §3.2 表 3、图 1 |
| `figures/` | 出版级图表（PNG） | 图 2、图 3 |

---

## 01_numeric_validation — 数值验证（48/48）

**文件：** `numeric_validation_raw.tsv`（48 数据行）

| 列 | 说明 |
|---|---|
| `pose` | 刚体变换：identity / translation / rotation / composite |
| `method` | 评分函数名（12 个可对比函数） |
| `system` | 测试体系：2oob（蛋白）/ 1azp（蛋白-DNA） |
| `python_score` | Python LightDock 0.9.4 参考分值 |
| `rust_score` | LKlight 分值 |
| `diff_abs` | \|python − rust\| 绝对差（容差 10⁻⁶；实测最大 4.85×10⁻⁷） |
| `status` | OK = 通过 |

**辅助：** `numeric_validation_v2.py`（复现脚本）

---

## 02_bm5_30case — 30 例 BM5 CAPRI 成功率

**文件：**
- `summary.tsv` — 逐复合物 top-N 成功率（论文表 4 数据源）
- `cumulative.tsv` — 逐检查点累计视图（同字段）

**字段（两文件一致）：**

| 列 | 说明 |
|---|---|
| `case` | BM5 复合物 PDB 代码（30 例） |
| `engine` | `lk`（LKlight）/ `py`（Python LightDock 0.9.4） |
| `n_solutions` | 有效模型数（top-N 分母） |
| `top1 … top100` | top-N 集内含 ≥1 个 CAPRI 可接受模型（Fnat≥0.1 且 L-RMSD≤10 Å）的成功率 |
| `status` | `ok` = 运行完成 |

**关键统计（论文 §3.3）：** 平均 \|Δ\| ≤ 0.020（全区间）；top-100 命中 lk 15/30 vs py 11/30（both 8，lk-only 7，py-only 3）；McNemar p=0.34；Wilcoxon p=0.12；并集成功率 18/30。

**提取示例：**
```python
import pandas as pd
df = pd.read_csv("02_bm5_30case/summary.tsv", sep="\t")
lk = df[df.engine == "lk"]; py = df[df.engine == "py"]
print(lk[["top5","top100"]].mean())   # 平均成功率
```

---

## 03_system_tests — 系统测试（档位 × 评分 × 场景）

**文件：**
- `tiers_eval.tsv` / `tiers_summary.tsv` — 档位 T1–T5（2X9A, fastdfire）
- `scoring_eval.tsv` / `scoring_summary.tsv` — 10 评分函数（2X9A）
- `scenes_eval.tsv` / `scenes_summary.tsv` — 5 应用场景
- `report_data.json` — 图表数据源（含耗时、原子数）

**统一字段：**

| 列 | 说明 |
|---|---|
| `tag` | 实验标识（tier_T1…/scoring_xxx/scene_xxx） |
| `n_models` | 有效模型数 |
| `best_scoring` / `median_scoring` | 最优/中位评分值 |
| `top1 … top100` | CAPRI 成功率 |

**tiers 档位参数（见 `tiers_summary.tsv` / report_data.json）：**
T1 10×50×100、T2 10×100×100、T3 20×100×100、T4 20×200×100、T5 25×200×300（swarms×glowworms×steps）

**关键结论（论文 §3.4）：** ① 五档成功率完全一致（top-5=0.60）——固定表面点下参数扩展零增益；② 评分三梯队（fastdfire/dfire 0.60 > dfire2 0.40 > 其余偶发；vdw/cpydock 零命中）；③ 耗时随原子数超线性（O(N²)）。

---

## 04_win_crossplat — Windows 跨平台验证

**文件：**
- `REPORT.md` — 完整验证报告（110/110 检查点字节级一致，最大 \|Δ\|=0.0）
- `score_timings.csv` / `score_timings_ps.csv` — Windows 端性能计时（2 核 Xeon 限制下仅供参考）
- `numeric_win.csv` — Windows 端数值验证记录

**验证设计：** 2X9A 体系，官方表面点初始位置，10 swarm × 50 glowworm × 100 步 fastdfire；macOS arm64 vs Windows x86-64（交叉编译）共享同一输入，GSO 轨迹逐位一致。

---

## 05_performance — 性能基准

**文件：**
- `raw_timings_20260827_110033.tsv` — 原始计时（3 次均值 × 4 场景 × 3 引擎）
- `table2_full.tsv` — 早期多评分函数 × 复合物成功率表（48 例体系扩展视图）

**raw_timings 字段：** `scenario`（1PPE-pydock / 1PPE-dfire / 1AZP-dnaANM / 1PPE-cpydock）、`engine`（python / rust-orig / lklight）、`wallclock_ms`、`note`（rc=0 正常；rc=101 = Rust-orig DFIRE 崩溃）

**论文表 3（最终数据，ms；Mac mini M4 10 核，单 swarm、100 步、200 萤火虫、3 次均值；`raw_timings_20260828_1130.tsv`）：**

| 场景 | Python | Rust-orig | LKlight | LKlight/Py | LKlight/官方 |
|---|---|---|---|---|---|
| 1PPE dfire | 164,968 | 崩溃 | 1,535 | 107.5× | — |
| 1PPE vdw | 25,721 | 崩溃 | 1,370 | 18.8× | — |
| 1PPE cpydock | 30,313 | 7,567 | 2,143 | 14.2× | 3.5× |
| 1AZP dna+ANM | 50,607 | 崩溃 | 3,901 | 13.0× | — |

> **数据勘误（重要）**：早期 `raw_timings_20260827_110033.tsv` 与 `benchmark.sh` 的 Python 计时无效——
> 该脚本以 `|| true` 吞掉了 lightdock3.py 的报错（初始位置文件数不匹配、缺 `lightdock.scoring.pydock` 模块），
> 所谓 "Python 381–404 ms" 实为启动+报错耗时。有效计时须以 `lightdock3 -l 0` 单 swarm 方式重新测量
> （见 `raw_timings_20260828_1130.tsv`）。LightDock 0.9.4 无 `pydock` 评分模块，pyDock 家族以 `cpydock` 代表。

---

## 引用说明

论文中引用本数据集的表述（投稿时可直接使用）：

> All benchmark data and the full reproducibility pipeline are available in the
> repository (`verify_equivalence/`), including per-model CAPRI metrics, the
> Windows cross-platform bit-identity validation, and the system-robustness
> test results.
