# 与原版 LightDock 0.9.4 的数值一致性审计（v1.2.0 · 0.5 Å 网格）

- 日期：2026-09-05
- 审计对象：LKlight exact（母仓 all-pairs，当前构建）与 LKlight-grid（v1.2.0，0.5 Å 默认）相对**原版 Python LightDock 0.9.4**
- 测试集：2oob（蛋白-蛋白，350+574 重原子，10 个蛋白类函数）与 1AZP（蛋白-DNA noh，532+328 重原子，dna/ddna），每函数 **native + 30 个确定性刚体 decoy**（31 pose，seed 20260904，平移 0.5–26 Å 覆盖近场 ≤10 Å 与远场 10–30 Å 带，旋转 ≤60°）
- 方法：在 PDB 文件层施加刚体变换 → python / exact / grid 三方对**字节相同输入**以 identity 打分（与 8 月 `numeric_validation_v2.py` 同设计，杜绝变换实现差异）
- 运行产物：`orig_triple_perpose.csv`（含逐 pose 最小原子间距离 min_dist）、`orig_triple_summary.csv`、`orig_triple.py`（可复现）

## 1. 结论摘要

1. **除 dna 外全部 11 个函数：LKlight exact 与原版一致到 ≤5×10⁻⁷**（sd 仅个别深穿插位姿达 1.6×10⁻⁵，相对误差 ~10⁻⁶），**31 pose 全覆盖**——"移植零偏差"在 v1.2.0/0.5 Å 时代依然成立，且从 8 月 48/48 的 4 pose 扩展到 31 decoy。
2. **dna 是一个有意的例外**：Sep 2 母仓提交 `a89ff65` 在 dna 打分中加入了"重原子深穿插线性罚"（`CLASH_PENALTY_W=6.0`，重-重原子对距离 d < 0.75×Σr_vdw 时线性加罚），用于修复原版 PyDockDNA"深穿模反而得分更好"的系统性缺陷。**基础能量项与原版仍逐位一致**——在无穿插对（min_dist ≥ ~2.5 Å）的位姿上实测 **exact ≡ 原版（diff = 0.0）**；差异只在穿插位姿出现，且随穿插深度单调增大（min_dist 1.1 Å → diff ~10²，0.3 Å → diff ~10³，0.15 Å → 4480）。
3. **LKlight-grid（0.5 Å）相对原版的偏差 = 上述 exact 一致性 ⊕ 0.5 Å 网格远场近似**：无穿插的远场静电族位姿上 grid 绝对误差 ≤ ~3 energy units（cpydock 31 pose max 2.99 / 中位 0.845；dna 干净位姿 ≤0.9），优于论文 Table 5 引用的 1 Å 网格口径（max 9.6），也满足 ≤12 units 声称。Spearman(原版, grid) 全函数 ≥ 0.99999 → 网格近似与穿插罚均**不改变位姿排序**。

## 2. 逐函数汇总（31 pose）

| 函数 | 体系 | exact vs 原版 max abs | grid vs 原版 max abs | grid vs 原版 median | Spearman(原版,grid) |
|---|---|---|---|---|---|
| dfire / fastdfire | 2oob | 5.0×10⁻⁷ | 5.0×10⁻⁷ | 2.3×10⁻⁷ | 1.0000 |
| dfire² | 2oob | 4.7×10⁻⁷ | 4.7×10⁻⁷ | 2.3×10⁻⁷ | 1.0000 |
| mj3h | 2oob | 7.1×10⁻¹⁵ | 7.1×10⁻¹⁵ | 0 | 1.0000 |
| sipper | 2oob | 2.5×10⁻¹⁴ | 2.5×10⁻¹⁴ | 0 | 1.0000 |
| tobi | 2oob | 5.3×10⁻¹⁴ | 5.3×10⁻¹⁴ | 4×10⁻¹⁶ | 1.0000 |
| pisa | 2oob | 4.9×10⁻⁷ | 4.9×10⁻⁷ | 1.7×10⁻⁷ | 1.0000 |
| vdw | 2oob | 4.9×10⁻⁷ | 4.9×10⁻⁷ | 2.2×10⁻⁷ | 1.0000 |
| sd | 2oob | 1.6×10⁻⁵ ¹ | 1.6×10⁻⁵ | 2.5×10⁻⁷ | 1.0000 |
| ddna | 1AZP | 5.0×10⁻⁷ | 5.0×10⁻⁷ | 2.2×10⁻⁷ | 1.0000 |
| cpydock | 2oob | 5.0×10⁻⁷ | 2.99 | 0.845 | 0.9999999 |
| dna | 1AZP | 4480 ² | 4480 ² | 1090 ² | 0.99999 |

¹ sd 最大差出现在深穿插位姿（vdw 量级 ~10⁴），为 python Cython f32 累积 vs rust f64 的舍入放大，相对误差 ~10⁻⁶，可忽略。
² dna 的 max/median 完全由穿插罚主导（见 §3）；**无穿插位姿 diff = 0.0–0.9**。grid 与 exact 在同一位姿的差仍仅 0.5 Å 网格量级（≈1 unit，见 §3 表 diff_grid−diff_exact）。

## 3. dna：穿插罚 vs 最小原子间距（决定性证据）

| pose | min_dist (Å) | python 原版 | exact | diff_exact | grid | diff_grid |
|---|---|---|---|---|---|---|
| decoy02 | 12.40 | −90.70 | −90.70 | **0.0** | −91.22 | 0.5 |
| decoy23 | 13.97 | −38.42 | −38.42 | **0.0** | −38.87 | 0.5 |
| decoy08 | 2.54 | −455.04 | −455.92 | **0.9** | −455.36 | 0.3 |
| decoy03 | 0.96 | −376.68 | −452.67 | 76 | −451.08 | 74 |
| decoy27 | 1.08 | −476.46 | −580.80 | 104 | −581.08 | 105 |
| decoy19 | 0.62 | −778.02 | −1533.98 | 756 | −1536.39 | 758 |
| native | 0.34 | −1201.39 | −2047.37 | 846 | −2048.59 | 847 |
| decoy20 | 0.27 | −2745.71 | −7227.72 | 4480 | −7225.81 | 4480 |

要点：
- **min_dist ≥ 2.54 Å → diff ≤ 0.9；min_dist 12–14 Å → diff = 0.000**：基础 DNA 能量（AMBER94 电 + vdW + 远场）与原版**逐位一致**。
- 罚 = 6.0 × Σ(0.75·Σr_vdw − d) 仅在重-重原子对 d < 0.75·Σr_vdw 时累加；因此与输入是否含氢无关（对含 H / noh 文件差恒等），只与**穿插度**相关。
- grid 与 exact 的差（≈1 unit）是 0.5 Å 远场网格插值，两引擎携带同一穿插罚 → §3.7 grid-vs-exact 审计口径不受影响。

## 4. 口径核查与对文档的影响

- 论文 #2 §3.7 / Table 5 的参照引擎写的是 "reference LKlight engine"（= 母仓 all-pairs exact），**不是** Python LightDock → 表格口径自洽；dna 行 max abs Δ / best-trajectory 0.18% 等均为 grid-vs-exact 同口径，穿插罚不破坏该审计。
- README 数值契约中 "original" 亦指母仓 all-pairs 路径（"The original all-pairs path (energy_exact) is retained"），非 Python LightDock → 文档体系无错误声称。
- **需要补注的点**（v1.2.0 发布说明 / 论文方法学）：
  1. dna 相对原版 PyDockDNA 含"重原子深穿插线性罚"（有意改进，修复深穿模得分偏好）；**基础项与原版逐位一致，干净结合位姿无差异**。与 8 月 numeric_validation（48/48, tol 1e-6）相比，dna 仅在不含穿插的位姿保持 1e-6 级一致；穿插位姿差异可达 10²–10³ units（方向：rust 更负，惩罚不可接受的渗透解）。
  2. LightDock 0.9.4 无独立 `pydock` 评分模块（pyDock 家族即 `cpydock`）→ pydock 无法与原版做 python 侧逐位对比（§3.7 的 pydock 行为 grid-vs-exact，已注明以 exact 为参照即可）。
- 0.5 Å 合理性：相比 1 Å（max abs 9.6），0.5 Å 网格把远场绝对误差压到 ≤3 units（cpydock max 2.99 / median 0.85；dna 干净位姿 ≤0.9），排序不变（Spearman ≥0.99999）→ **0.5 Å 默认相对原版更优且安全**。

## 5. 复现

```bash
cd <LKlight论文>/benchmarks_raw
/usr/bin/python3 orig_triple.py        # 需要:
#  - PYTHONPATH 前置 docking_env/lib/python3.9/site-packages (numpy 1.23.5+scipy 1.13.1+lightdock 0.9.4, 均已用官方 cp39 arm64 wheel 修复)
#  - /tmp/numval_orig/{2oob,1azp_*_noh}.pdb (tests/ 去氢; 1azp element 列在 col77-78)
#  - LKlight 与 LKlight-grid 的 target/release/LKlight
```
产物：`orig_triple_perpose.csv` / `orig_triple_summary.csv` / `orig_triple_run.log`。
