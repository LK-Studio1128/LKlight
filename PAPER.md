# LKlight: A High-Performance Rust Reimplementation of the LightDock Glowworm Swarm Optimization Docking Engine

**LKlight：基于 Rust 语言的 LightDock 萤火虫群优化对接引擎高性能再实现**

> **Version:** 1.0.0 | **Base:** LightDock 0.9.4 (Python + lightdock-rust) | **Repository:** https://github.com/LK-Studio1128/LKlight | **License:** GPL-3.0-or-later
> **Binary:** `LKlight` | **Key features:** 12 scoring-function families / 13 CLI method names · ANM support · rayon parallel outer loop · SIMD-friendly hot paths · thread-local scratch reuse · BufWriter I/O · macOS arm64 / Linux x86-64 / Windows x86-64

---

## Abstract

Molecular docking is a cornerstone of structure-based drug design and protein–protein interaction (PPI) analysis. LightDock [1,2], developed at the Barcelona Supercomputing Center, is an open-source docking framework built upon the Glowworm Swarm Optimization (GSO) meta-heuristic [3]. It provides 12 scoring functions spanning statistical potentials (DFIRE [4], DFIRE2 [5]), physicochemical potentials (PyDock [6], cpyDOCK), and biophysical potentials (PISA [8], SIPPER [9], MJ3h [10], TOBI [11]), and natively supports backbone flexibility through Anisotropic Network Models (ANM) [7]. While LightDock's Python implementation offers algorithmic breadth, interpreter overhead limits its throughput in large-scale virtual screening contexts.

We present **LKlight v1.0**, a complete reimplementation of the LightDock computational core in safe Rust [12]. LKlight first corrects four critical defects in the prior Rust baseline (`lightdock-rust`): (i) runtime panics from missing DFIRE parameter files, eliminated by binary embedding; (ii) erroneous ANM stride computations across multiple scoring modules causing out-of-bounds reads; (iii) unguarded panics on non-standard residues in DFIRE; and (iv) a mismatched atom-count assertion disabling all ANM-enabled runs. Building on this corrected foundation, LKlight applies a multi-tier optimization strategy comprising `rayon`-based [13] parallelization of the receptor-atom outer loop, SIMD-friendly contiguous hot loops compatible with portable release baselines, `thread_local!` scratch-buffer reuse eliminating per-step heap allocations, a 3D spatial hash grid for the SD scoring function (9 Å cutoff, genuine O(N²)→O(N) sparsification), and `BufWriter` I/O batching. We also provide a rigorous quantitative analysis of why analogous grids regress performance for scoring functions with large cutoff radii (≥15 Å).

Benchmarks on macOS arm64 (Apple Silicon, 200 glowworms, 100 steps, *n* = 3 replicates) demonstrate speedups of **3.0–25.5× over Python** and **26.5–307× over the prior Rust baseline** across four representative docking scenarios. Correctness is guarded by the public `cargo test --lib` suite and by development-stage numerical integration comparisons against the Python reference implementation. LKlight is released as a GPL-3.0-or-later derivative of LightDock, with full source code availability on GitHub and pre-built binaries distributed separately as release assets.

---

**中文摘要**

分子对接是基于结构的药物设计与蛋白质-蛋白质相互作用分析的核心计算工具。LightDock 是一个基于萤火虫群优化（Glowworm Swarm Optimization, GSO）算法的开源分子对接框架，提供 12 种评分函数，支持蛋白质-蛋白质、蛋白质-DNA 及抗体-抗原对接，并原生支持各向异性网络模型（ANM）柔性。然而，原始 Python 实现受限于解释型语言的内在性能上限，在大规模虚拟筛选场景下计算开销显著。

本文介绍 **LKlight v1.0**，一个将 LightDock 核心引擎从 Python/NumPy 完整迁移至 Rust 的实现。主要贡献包括：

1. 将 12 类评分函数（13 个命令行方法名，其中 `fastdfire` 为 `dfire` 的兼容别名）完整移植至 Rust，修复了原 Rust 基线版本中存在的多项 Bug（DFIRE 参数文件缺失崩溃、ANM stride 计算错误等）；
2. 实施多层次性能优化：预计算 `sqrt_vdw_charges`（消除热路径 sqrt 调用）、`thread-local` 坐标/界面 Vec 复用（消除每步堆分配）、GSO 运动阶段并行化（`rayon par_iter_mut`）、空间网格剪枝（sd.rs，O(N²)→O(N)）、BufWriter 减少文件写 syscall；
3. 通过对四类测试场景（1PPE pydock、1PPE dfire、1AZP DNA+ANM、1PPE cpydock）的系统基准测试，定量揭示各优化项的实际效果，并分析空间网格优化在大截断距离场景下出现性能回退的根因。

基准测试（macOS arm64，swarm_0，200 glowworms，100 步，3次均值）表明：

- **pydock**（蛋白质-蛋白质，30Å 截断）：LKlight **290ms** vs Python **858ms**，**快 3.0×**；vs Rust-orig **7693ms**，**快 26.5×**
- **dna+ANM**（蛋白质-DNA，ANM 模式）：LKlight **46ms** vs Python **760ms**，**快 16.5×**；vs Rust-orig **14142ms**，**快 307×**
- **cpydock**（含去溶剂化）：LKlight **44ms** vs Python **844ms**，**快 19.2×**；vs Rust-orig **7158ms**，**快 163×**
- **dfire**（哈希表统计势）：LKlight **33ms** vs Python **840ms**，**快 25.5×**（Rust-orig 因缺少外部参数文件而**运行崩溃**，约 6ms 为启动崩溃时间）

性能突破的关键是多项联合优化：(1) 撤回 HashMap 空间网格（G3/G4 回退）；(2) 用 `rayon par_iter` 将受体原子外循环并行化（H1，pydock/dna/cpydock）；(3) 将同一并行化模式应用于 dfire/dfire2/sd（I1-I3）；(4) 将热路径重构为 SIMD 友好的连续数组和简单内循环，同时发布配置采用便携 CPU baseline，benchmark 构建可按需启用 native 优化。产出的单文件二进制 `LKlight` 支持 macOS arm64、Linux x86-64 和 Windows x86-64。

---

## 1. Introduction

### 1.1 分子对接背景

分子对接是计算结构生物学的核心工具，广泛应用于蛋白质-蛋白质相互作用预测、蛋白质-DNA 对接以及基于结构的药物设计。当前主流对接工具包括 HADDOCK [15]（数据驱动共测定）、AutoDock Vina [16]（小分子对接）以及 LightDock [1,2]（萤火虫群优化，多评分函数），各自在特定场景下具有独特优势。

LightDock 由 Jiménez-García 等人 [1] 在巴塞罗那超级计算中心开发，其核心算法为萤火虫群优化（GSO）[3]。与传统 Monte Carlo 或遗传算法不同，GSO 模拟萤火虫种群通过相互感知萤光素浓度（luciferin）决策移动方向，能够并发维持多个局部极值，理论上更适合高维构象空间的多模态搜索。 LightDock v2.0 [2] 进一步引入信息驱动约束，支持将实验数据（除NMR化学位移等）作为对接引导。

LightDock 的关键特性包括：

| 特性 | 说明 |
|------|------|
| **GSO 引擎** | 每个 swarm 运行 N_g 个 glowworm，步进 T 步，luciferin 更新 + 概率移动选择 |
| **多评分函数** | 12 种覆盖物理势（DFIRE、VDW）、知识势（SIPPER、MJ3h）、物理化学势（PyDock、SD）和生物物理势（PISA、TOBI、dDNA） |
| **ANM 支持** | 通过各向异性网络模型（ANM）描述受体/配体骨架柔性，模式向量 nmodes 以 1D 数组内嵌于 GSO 状态向量 |
| **多种对接类型** | 蛋白质-蛋白质、蛋白质-DNA、抗体-抗原、跨膜蛋白-膜蛋白 |

原始 Python 实现通过 NumPy 广播运算实现配对能量计算的隐式 SIMD 向量化，在小规模 swarm 数量下具有合理性能，但解释型开销在每步 GSO 框架（luciferin 更新、邻居搜索、运动阶段）中仍不可忽视。

### 1.2 已有 Rust 基线版本的问题

在本工作开展之前，已存在一个 LightDock 的 Rust 基线版本（`lightdock-rust`，包含于 `lightdock-macos-arm64` 发行包）。通过系统测试，我们发现该基线版本存在以下问题：

1. **DFIRE 参数文件缺失**：基线版本在运行 `dfire` 评分时会搜索外部参数文件，文件不存在时直接 `panic`，导致 DFIRE/DFIRE2/DDNA 评分函数完全无法使用；
2. **ANM stride 计算错误**：`pisa.rs`、`ddna.rs`、`cpydock.rs`、`sd.rs` 等模块在处理 ANM 模式向量时，stride 计算不当，导致数组越界或错误计算；
3. **未知残基 panic**：`dfire.rs` 在遇到标准 20 种氨基酸以外的残基（如配体、修饰残基）时直接崩溃，而非优雅降级；
4. **不必要的 atom_count() 断言**：`simulator.rs` 中的断言与 ANM 模式向量长度不匹配，导致含 ANM 的对接任务失败。

### 1.3 本文贡献

| 贡献类别 | 具体内容 |
|----------|---------|
| **Bug 修复** | DFIRE 参数嵌入（消除外部文件依赖）；ANM stride 统一为 `nmodes.len()/(3×n_modes)`；未知残基返回 999 而非 panic；移除不兼容的 atom_count() 断言 |
| **评分函数完整性** | 12 类评分函数（13 个命令行方法名，`fastdfire` 为 `dfire` 兼容别名）均有 Rust 实现；公开仓库保留 `cargo test --lib` 单元测试与轻量 PDB 夹具，开发阶段另以 Python 参考实现进行综合数值对比 |
| **内存优化** | `thread_local!` 坐标/界面向量复用；GSO 运动阶段 pos_scratch/rot_scratch 字段复用；`qt::rotate()` 返回 `[f64;3]` 消除堆分配 |
| **计算优化** | `sqrt_vdw_charges` 预计算；sd.rs 9Å 空间网格 O(N²)→O(N)；GSO 运动阶段 rayon 并行化 |
| **I/O 优化** | `swarm.rs save()` 使用 `BufWriter` 批量写出减少 syscall |
| **空间网格分析** | 定量分析 10Å/±3格方案在大截断距离（30Å）下出现性能回退的根因，为后续优化提供指导 |

---

## 2. Theory

### 2.1 萤火虫群优化（GSO）算法

GSO 的核心迭代包含三个阶段：

**萤光素更新（Luciferin update）**

$$\ell_i(t+1) = (1-\rho)\,\ell_i(t) + \gamma\, J(x_i(t))$$

其中 $\rho$ 为萤光素衰减率，$\gamma$ 为增强系数，$J(x_i)$ 为当前位置的评分函数值。

**邻居感知与概率选择（Neighbor selection）**

$$N_i(t) = \{j \mid \|x_j - x_i\| < r_i^d(t),\ \ell_j(t) > \ell_i(t)\}$$

$$P(i \to j) = \frac{\ell_j(t) - \ell_i(t)}{\sum_{k \in N_i} \ell_k(t) - \ell_i(t)}$$

**运动更新（Movement update）**

$$x_i(t+1) = x_i(t) + s\,\frac{x_j - x_i}{\|x_j - x_i\|}$$

在 6-DOF（3 平移 + 3 旋转）分子对接扩展中，旋转分量使用四元数球面线性插值（SLERP）。

**视野范围自适应（Vision range update）**

$$r_i^d(t+1) = \min\!\left(r_s,\, \max\!\left(0,\, r_i^d(t) + \beta(n_t - |N_i(t)|)\right)\right)$$

### 2.2 评分函数

本实现支持 12 类评分函数、13 个命令行方法名（`fastdfire` 为 `dfire` 的兼容别名），涵盖三类势能：

**统计势（Statistical Potentials）**

DFIRE 基于原子对距离频率的参考态模型：
$$E_{\text{DFIRE}} = \sum_{i<j} \Delta E_{\text{DFIRE}}(\Delta_r^{ij},\ t_i,\ t_j)$$

**物理化学势（Physicochemical）**

PyDock 包含静电和 van der Waals 两项：
$$E_{\text{pyDOCK}} = -\left(\frac{F}{\varepsilon}\sum_{i,j}\frac{q_i q_j}{r_{ij}^2}\cdot\mathbf{1}[r_{ij} \leq r_{\text{elec}}] + \sum_{i,j}\varepsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - 2\left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right]\cdot\mathbf{1}[r_{ij} \leq r_{\text{vdw}}]\right)$$

**结构去溶剂化（Desolvation）**

cpyDOCK 在 PyDock 基础上引入原子溶剂可及面积（ASA）加权的去溶剂化项：
$$E_{\text{solv}} = \sum_{i}\text{ASA}_i \cdot d_i \cdot \min\!\left(-10 d_{\min,i} + 65,\ \text{ASA}_i\right) \cdot \mathbf{1}\!\left[d_{\min,i} \leq d_{\text{solv}}\right]$$

### 2.3 ANM 柔性处理

各向异性网络模型（ANM）以刚体运动叠加正则模式的方式引入主链柔性：

$$\mathbf{x}_{\text{anm}}(k) = \mathbf{x}_0(k) + \sum_{m=1}^{M} q_m \cdot \mathbf{v}_m(k),\quad k=1,\ldots,N_{\text{anm}}$$

其中 $\mathbf{v}_m(k)$ 为第 $m$ 个正则模式在原子 $k$ 处的方向向量（$3N_{\text{anm}}$ 维），$q_m$ 为 GSO 状态向量中对应的模式振幅分量。模式向量以 stride = $N_{\text{anm}} \times 3$ 存储为 1D `Vec<f64>`，索引为：

$$b = m \cdot N_{\text{anm}} \cdot 3 + k \cdot 3, \quad \text{stride} = \frac{\text{nmodes.len()}}{3 \cdot n\_modes}$$

这一正确的 stride 计算是修复 ANM Bug 的关键（见第 3.2 节）。

---

## 3. Methods

### 3.1 完整评分函数移植与命令行方法名

12 类评分函数均经过独立 Rust 实现，并在开发阶段与 Python 参考值进行数值验证；公开仓库中保留可由 `cargo test --lib` 运行的核心单元测试，以及 `tests/` 下的轻量 PDB 夹具用于示例和烟雾验证：

| 评分函数 | 类别 | 截断距离 | ANM 支持 | 关键特性 |
|---------|------|---------|---------|---------|
| `dfire` / `fastdfire` | 统计势 | 15 Å | ✓ | DFIRE 参数内嵌，残基对距离索引 |
| `dfire2` | 统计势 | 15 Å | ✓ | DFIRE2 参数内嵌 |
| `dna` | 物理化学 | 30 Å (elec) / 10 Å (vdw) | ✓ | DNA-蛋白质专用参数 |
| `mj3h` | 知识势 | 7.5 Å | ✓ | MJ 残基接触矩阵 |
| `pydock` | 物理化学 | 30 Å (elec) / 10 Å (vdw) | ✓ | Coulomb 静电 + LJ vdW |
| `cpydock` | 物理化学+溶剂 | 30 Å (elec) / 6.4 Å (solv) | ✓ | pydock + 去溶剂化 |
| `sd` | 溶剂化 | 9 Å | ✓ | ASA 加权接触去溶剂化 |
| `vdw` | 物理 | 10 Å | ✓ | 纯 LJ van der Waals |
| `pisa` | 统计势 | 8.5 Å | ✓ | PISA 原子接触统计 |
| `sipper` | 统计势 | 8.5 Å | ✗ | 残基对接触势 |
| `tobi` | 统计势 | 12 Å | ✓ | TOBI 原子对势，消除 sqrt |
| `ddna` | 统计势 | 15 Å | ✓ | dDNA 统计势 |

其中 `parse_method()` 接受的命令行方法名为：

```
dfire fastdfire dfire2 dna mj3h pydock cpydock sd vdw pisa sipper tobi ddna
```

### 3.2 Bug 修复详情

#### Fix 1 — DFIRE 参数内嵌（消除崩溃）

**问题**：Rust 基线版本在 `dfire.rs` 中通过文件路径加载 DFIRE 参数矩阵，路径硬编码为相对于运行目录的固定位置。当从任意目录运行或部署为单一可执行文件时，参数文件不存在导致 `unwrap()` 崩溃：

```
thread 'main' panicked at src/dfire.rs:247:14:
Unable to open DFIRE parameters: Os { code: 2, kind: NotFound, message: "No such file or directory" }
```

**修复**：将 DFIRE、DFIRE2、DDNA 参数矩阵以 `include_bytes!` 或常量数组形式嵌入二进制，`new()` 函数直接从嵌入数据初始化，无需任何外部文件依赖。

#### Fix 2 — ANM Stride 统一

**问题**：多个评分函数（`pisa.rs`、`ddna.rs`、`cpydock.rs`、`sd.rs`）中，ANM 模式向量的 atom stride 计算为固定值（如 `3 * n_atoms`），与实际 nmodes 存储布局不符，导致 stride 与 `nmodes.len()/(3*n_modes)` 不匹配时出现越界访问或错误坐标。

**修复**：统一所有评分函数的 ANM stride 计算为：

```rust
let nm_n = if num_anm > 0 {
    nmodes.len() / (3 * num_anm)
} else { n_atoms };
```

同时加入 bounds guard `if i_atom >= nm_n { break; }` 防止越界。

#### Fix 3 — DFIRE 未知残基优雅降级

**问题**：`dfire.rs` 在遇到标准 20 种氨基酸以外的残基（非标准残基、小分子、修饰氨基酸）时，未能正确查表并直接 `panic`。

**修复**：`new()` 函数遇到未知残基时发出警告并 `continue` 跳过，而非 panic；`energy()` 对未知残基类型返回 999（最大惩罚值）而非崩溃。

#### Fix 4 — simulator.rs atom_count() 断言

**问题**：`simulator.rs` 中的断言 `assert_eq!(model.atom_count(), anm_atoms)` 与 ANM 只覆盖 Cα 原子（而非全原子）的事实不符，导致所有含 ANM 的对接任务失败。

**修复**：移除该断言，改用 stride 推导实际 ANM 原子数。

### 3.3 内存与计算优化（Session 3–6）

#### F1 — sqrt_vdw_charges 预计算

pydock、cpydock、sd、dna 的 VDW 能量热路径中包含 `(rec.vdw_charges[i] * lig.vdw_charges[j]).sqrt()`，每次调用均产生一次 `sqrt` 计算（~20 cycles on x86）。预计算 `sqrt_vdw_charges` 字段，在 `new()` 时一次性完成：

```rust
pub sqrt_vdw_charges: Vec<f64>,
// ...
sqrt_vdw_charges.push(vdw_charge.sqrt());
```

热路径替换为乘法：`rec.sqrt_vdw_charges[i] * lig.sqrt_vdw_charges[j]`。

#### F2 — sd.rs 9Å 空间网格（O(N²)→O(N)）

SD 评分函数的接触截断距离为 9 Å，远小于蛋白质尺寸，大量受体-配体原子对均超出截断。使用 `thread_local!` `HashMap<(i32,i32,i32), Vec<usize>>` 对配体坐标建立 3D 网格（cell = 9 Å），外循环只查询受体原子邻域 27 个格点：

$$\text{cell}(x) = \lfloor x / \text{CELL} \rfloor, \quad \text{search: } \Delta\in\{-1,0,+1\}^3 = 27\ \text{cells}$$

格点在每次能量评估前 `clear()` 并重建，通过 `or_default()` 复用已分配的内部 Vec，避免每步重新分配。

#### F3 — BufWriter 减少文件 I/O syscall

`swarm.rs save()` 函数将每步 gso_*.out 文件改为 `BufWriter<File>` 写出，将每行一次 `write` syscall 合并为批量写出：

```rust
let mut writer = BufWriter::new(File::create(path)?);
```

#### F4 — qt::rotate() 返回 [f64;3]

四元数旋转函数 `rotate()` 原返回 `Vec<f64>`（堆分配），改为返回 `[f64;3]`（栈上固定数组），消除热路径中的每次向量分配：

```rust
pub fn rotate(&self, v: [f64; 3]) -> [f64; 3] { ... }
```

#### F5 — glowworm.rs 栈上坐标

`Glowworm` 的 `translation` 字段从 `Vec<f64>` 改为 `[f64; 3]`，`move_towards` 内部的 delta_x 中间计算从 Vec 改为栈数组，消除 GSO 运动阶段的频繁小向量分配。

#### G1 — Swarm pos_scratch/rot_scratch 字段复用

`Swarm` 结构体新增 `pos_scratch: Vec<[f64; 3]>` 和 `rot_scratch: Vec<Quaternion>` 字段，在 `movement_phase()` 中复用（`resize + fill`），替代每步 `Vec::new() + push()` 模式：

```rust
pub struct Swarm<'a> {
    pub glowworms: Vec<Glowworm<'a>>,
    pos_scratch:   Vec<[f64; 3]>,
    rot_scratch:   Vec<Quaternion>,
}
```

同时将 `neighbors.into_iter()` 替代 `neighbors[i].clone()`，消除 $N_g \times k$ 次邻居列表克隆。

#### G2 — movement_phase 并行化

使用 `rayon` 将 GSO 运动阶段并行化：预生成 `Vec<f64>` 随机数序列（解决线程安全问题），通过字段级分借（`gws`、`pos_s`、`rot_s` 来自不同字段）实现 `par_iter_mut`：

```rust
let randoms: Vec<f64> = (0..n).map(|_| rng.gen()).collect();
let (gws, pos_s, rot_s) = (&mut self.glowworms, &self.pos_scratch, &self.rot_scratch);
gws.par_iter_mut()
   .zip(randoms.par_iter())
   .for_each(|(gw, &r)| {
       let nid = gw.select_random_neighbor(r) as usize;
       gw.move_towards(nid as u32, &pos_s[nid], &rot_s[nid], ...);
       gw.update_vision_range();
   });
```

#### G3/G4 → H1 — pydock/dna/cpydock 外循环并行化（最终方案）

初始尝试（G3/G4）为 pydock/dna/cpydock 引入 10Å/±3格空间网格，结果产生性能回退（HashMap 查询开销 >> 计算节省，见第 4.3 节）。

**H1（最终实现）：撤回 HashMap 网格，改用 rayon 并行化受体原子外循环。** 实现分两阶段：

```rust
// Phase 1: parallel ELEC + VDW (receptor atoms × ligand atoms)
let (total_elec_raw, total_vdw) = receptor_coords.par_iter().enumerate()
    .map(|(i, ra)| {
        let (mut ei, mut vi) = (0.0f64, 0.0f64);
        for (j, la) in lig_slice.iter().enumerate() {
            let d2 = dist2(ra, la);
            if d2 <= ELEC_DIST_CUTOFF2 {
                ei += (rec_ele[i] * lig_ele[j] / d2).clamp(ELEC_MIN_CUTOFF, ELEC_MAX_CUTOFF);
            }
            if d2 <= VDW_DIST_CUTOFF2 {
                vi += (rec_svdw[i] * lig_svdw[j] * lj6_kernel(rec_vdwr[i]+lig_vdwr[j], d2)).min(VDW_CUTOFF);
            }
        }
        (ei, vi)
    })
    .reduce(|| (0.0, 0.0), |(e1,v1),(e2,v2)| (e1+e2, v1+v2));

// Phase 2: sequential interface flags (INTERFACE_CUTOFF=3.9Å, fast)
for (i, ra) in receptor_coords.iter().enumerate() {
    for (j, la) in lig_slice.iter().enumerate() {
        if dist2(ra, la) <= INTERFACE_CUTOFF2 { iface_r[i]=1; iface_l[j]=1; }
    }
}
```

Phase 2（interface flags）占总时间 <5%（截断仅 3.9Å，大多数对跳过），并行 Phase 1 贡献全部性能提升。

#### H2 — SIMD 友好热路径与便携发布 baseline

```toml
[target.x86_64-unknown-linux-gnu]
rustflags = ["-C", "target-cpu=x86-64"]

[target.x86_64-pc-windows-msvc]
rustflags = ["-C", "target-cpu=x86-64", "-C", "target-feature=+crt-static"]
```

LKlight 的优化策略分为两层：源码层面将热路径改写为连续数组访问、简单内循环和更少分支，使 LLVM 更容易进行自动向量化；发布层面则采用 `target-cpu=x86-64` 等便携 baseline，保证 Linux/Windows 二进制能在更广泛机器上运行。对于本机 benchmark 或内部性能测试，可临时使用 `RUSTFLAGS="-C target-cpu=native"` 构建以释放 AVX2/FMA 或 NEON/ASIMD 等平台特性，但这不是公开 Release 二进制的默认配置。

---

## 4. Results

### 4.1 正确性验证

公开仓库中的可复现正确性检查包括 `cargo test --lib` 单元测试和 `tests/` 下轻量 PDB 夹具；开发阶段另使用 Python LightDock 参考实现完成全评分函数数值对比：

| 测试类别 | 测试数量 | 通过数量 |
|---------|---------|---------|
| `cargo test --lib`（公开单元测试） | 29 | **29 / 29** |
| `tests/` 轻量 PDB 夹具 | 4 files | **Present** |
| 开发阶段综合数值对比（全评分函数） | 160 | **160 / 160** |

浮点精度：G3/G4 空间网格引入的 FP 累加顺序变化导致数值差异 ~2×10⁻¹³，相对误差 < 10⁻¹²，在科学计算精度范围内完全可接受。单元测试断言已从精确相等更新为带容差比较（ε = 10⁻⁸）。

### 4.2 性能基准测试

**测试环境：** macOS arm64（Apple Silicon），单 swarm，200 glowworms，100 步，每项3次重复取均值

**测试场景：**
- **1PPE pydock**：胰蛋白酶-BPTI 复合物（1615 受体原子 × 221 配体原子 = 357K 对），PyDock 评分（ELEC 截断 30Å + VDW 截断 10Å），不含 ANM
- **1PPE dfire**：同上结构，DFIRE 统计势
- **1AZP dna+ANM**：转录因子-DNA 复合物（ANM 模式开启），DNA 评分函数
- **1PPE cpydock**：胰蛋白酶-BPTI，cpyDOCK 评分（含去溶剂化）

| 测试场景 | Python (ms) | Rust-orig (ms) | LKlight (ms) | LKlight/Py× | LKlight/Orig× |
|---------|------------|--------------|--------------|--------|----------|
| 1PPE pydock | 858 | 7,693 | **290** | **3.0×** | **26.5×** |
| 1PPE dfire | 840 | **CRASH** ¹ | **33** | **25.5×** | N/A |
| 1AZP dna+ANM | 760 | 14,142 | **46** | **16.5×** | **307×** |
| 1PPE cpydock | 844 | 7,158 | **44** | **19.2×** | **163×** |

> ¹ Rust-orig dfire 因外部参数文件缺失而在运行时 panic（exit code 101），~6ms 为启动+崩溃时间，非真实计算。

**平台：** macOS arm64（Apple Silicon），swarm_0，200 glowworms，100 步，3 次重复取均值。

**结论：**
- **pydock**：LKlight **3.0× 快于 Python**，**26.5× 快于 Rust-orig**
- **dna+ANM**：LKlight **16.5× 快于 Python**，**307× 快于 Rust-orig**
- **cpydock**：LKlight **19.2× 快于 Python**，**163× 快于 Rust-orig**
- **dfire**：LKlight **25.5× 快于 Python**（Rust-orig 崩溃，LKlight 为可用 Rust 实现，并且全面超越 Python）

### 4.3 G3/G4 空间网格回退分析与修复（H1）

**回退根本原因：** 10Å/±3格方案中，每个受体原子需执行 7³ = **343 次 HashMap 查询**（多数返回空）。对 1PPE 配体（221 原子，~14 个非空格点），实际有效查询比例约 14/343 = **4%**。

$$\text{每受体原子开销}_{\text{grid}} = 343 \times t_{\text{HashMap}} \approx 343 \times 50\text{ns} = 17\mu\text{s}$$
$$\text{每受体原子开销}_{\text{O(N²)}} = 221 \times t_{\text{array}} \approx 221 \times 2\text{ns} = 0.44\mu\text{s}$$

网格方案 HashMap 查询开销（~17μs）比直接数组遍历（~0.44μs）高约 **38×**。

**F2（sd.rs, 9Å grid）有效的原因：** SD 评分截断为 9Å，仅查询 **27 个格点**（3³），且大多数配体原子在截断之外，实际计算量大幅减少。该优化仍保留。

**H1 修复：** 撤回 G3/G4 中的 HashMap 网格，恢复简洁 O(N²) 内循环，同时添加 `rayon par_iter` **并行化受体原子外循环**（分两阶段：并行 ELEC+VDW，顺序 interface flags）。此修复将 pydock 从 Rust-orig 7693ms 降至 290ms（26.5× 提升）。

**I1/I2（dfire/dfire2 Session 8 新增）：** dfire 的原始 HashMap 空间网格（CELL=15Å，±1=27格）与 G3/G4 存在相同根因：15Å 网格覆盖 1PPE 整个受体，几乎不剪枝，但带来 27 次 HashMap 查询/配体原子固定开销。移除网格后改用 rayon 并行受体原子外循环，dfire 从 935ms 降至 **33ms**（35× 提升，超越 Python 25.5×）。

**I3（sd.rs Session 8 新增）：** sd.rs 的 9Å 空间网格（3³=27格）截断远小于蛋白质尺寸，真正稀疏有效。在受体原子外循环保留网格查询（Phase 1 并行），Phase 2 顺序更新 interface flags。

### 4.4 并行 + SIMD 联合优化分析（H2 + I1/I2/I3）

**H2：SIMD 友好热路径 + 可选 native benchmark 构建。** LKlight 的核心改动不是依赖不可移植的默认编译参数，而是把热路径改成适合 LLVM 自动向量化的形式：连续坐标数组、简单距离平方计算、较少临时分配和较少虚调用。公开发布二进制使用便携 CPU baseline；本机 benchmark 可使用 `target-cpu=native` 观察硬件上限。

$$\text{总加速} \approx N_{\text{cores}} \times W_{\text{SIMD}} = 8 \times 4 = 32\times \quad\text{（理论峰值）}$$
$$\text{实测（pydock）} = \frac{7693\text{ms}}{290\text{ms}} = 26.5\times \approx \text{理论值的 83\%}$$

**I1/I2 — dfire/dfire2 H1 化（移除 HashMap 网格）：** 将 dfire 的 HashMap 空间网格（CELL=15Å，±1=27格）替换为 rayon 并行受体原子外循环，全量 O(N²) 内循环。效果：dfire **35× 提升**（935ms → 33ms），超越 Python 25.5×。dfire 原始 HashMap 网格的失效原因与 G3/G4 相同：15Å 网格覆盖整个受体，几乎不剪枝，但带来 27 次 HashMap 查询的固定开销。

**I3 — sd.rs 并行化（保留 9Å 网格）：** sd.rs 的 9Å 网格仍然有效（3³=27格，短截断真正稀疏），在并行基础上保留网格进一步减少工作量。

**dfire DFIRE 表查询为何可 SIMD 化：** 内循环纯算术（距离计算 + 整数 bin 索引 + 数组访问），LLVM 可向量化距离计算部分，表查找部分为依赖 gather，在 arm64 NEON 上部分向量化。实测加速（25×）超出理论单路 SIMD 预期，说明 rayon 并行化是主要贡献（8 核 × ~3× SIMD = ~24×）。

---

## 5. Implementation

### 5.1 代码结构

```
src/
├── bin/
│   ├── lightdock.rs        # 统一入口（setup/run/rank/top/score/pipeline 等子命令）
│   ├── lightdock-rust.rs   # 兼容上游 run 入口
│   ├── lightdock-setup.rs  # setup 辅助入口
│   ├── lgd_rank.rs         # rank 辅助入口
│   └── lgd_generate_conformations.rs
├── swarm.rs            # Swarm 结构 + GSO 引擎（G1/G2 优化）
├── glowworm.rs         # Glowworm 结构 + 运动/概率（栈上坐标）
├── qt.rs               # 四元数 + SLERP（返回 [f64;3]）
├── simulator.rs        # ANM 应用 + 坐标变换
├── scoring.rs          # Score trait + satisfied_restraints + membrane_intersection
├── pydock.rs           # PyDock 评分（rayon 并行化 ✓，H1）
├── cpydock.rs          # cpyDOCK 评分（rayon 并行化 ✓，H1）
├── dna.rs              # DNA 评分（rayon 并行化 ✓，H1）
├── sd.rs               # SD 评分（F2 9Å 网格 ✓）
├── dfire.rs            # DFIRE/DFIRE2 评分（参数内嵌 ✓）
├── ddna.rs             # dDNA 评分
├── pisa.rs             # PISA 评分（空间索引 ✓）
├── tobi.rs             # TOBI 评分（空间索引 + sqrt 消除 ✓）
├── vdw.rs              # VDW 评分
├── mj3h.rs             # MJ3h 残基接触势
└── sipper.rs           # SIPPER 残基接触统计
```

代码总规模约 **8,100 行 Rust 源码**（不含 data/ 参数文件）。

### 5.2 统一 CLI 功能

公开发布的 `LKlight` 单二进制入口覆盖 LightDock 常用工作流、分析工具和辅助格式转换：

| 子命令 | 功能 |
|--------|------|
| `setup` | 从受体/配体 PDB 生成 `setup.json`、`initial_positions_*.dat`、swarm 初始目录；支持 `--anm` 与约束文件 |
| `run` | 对指定 swarm 初始位置运行 GSO 优化 |
| `generate` | 根据 GSO 输出生成 top 构象 PDB；支持 ANM-aware 坐标重构 |
| `cluster` | 对 GSO 输出进行 DBSCAN 式聚类 |
| `rank` / `rank_swarm` | 汇总全部 swarm 或逐 swarm 排名 |
| `top` | 从 ranking 文件生成 Top-N PDB |
| `filter` | 根据 restraints 文件过滤 ranking |
| `gso_to_csv` | 将 ranking / GSO 输出转换为 CSV |
| `move_anm` | 基于 ANM 模式生成柔性构象 |
| `score` | 对给定受体/配体 PDB 进行单点评分，可传入平移/四元数 |
| `diameter` | 计算 PDB 结构直径 |
| `trajectory` | 从某个 glowworm 的 GSO 轨迹生成逐步 PDB |
| `map_contacts` | 从对接构象映射受体-配体接触 |
| `reference_points` | 计算或保存结构参考点 |
| `pipeline` | 一条命令完成 setup、run、rank、top 的自动化流程 |

### 5.3 线程安全与并行策略

- **GSO 邻居搜索**：`par_iter()` 只读并行（`self.glowworms` 共享引用），无锁
- **GSO 运动阶段**：`par_iter_mut()` 修改每个 `Glowworm`，通过字段级分借（`glowworms` + `pos_scratch` + `rot_scratch`）满足借用检查
- **评分函数**：`thread_local! { static SCRATCH: RefCell<...> }` 确保每线程独立缓存，无竞争
- **随机数**：运动阶段随机数在并行前预生成（`StdRng`），保证确定性可复现

### 5.4 平台支持与二进制发布

| 平台 | 状态 |
|------|------|
| macOS arm64（Apple Silicon）| 完全支持；Release 资产建议命名 `LKlight-macos-arm64.tar.gz` |
| Linux x86-64 | 完全支持；静态/便携二进制 Release 资产建议命名 `LKlight-linux-x86_64.tar.gz` |
| Windows x86-64 | 完全支持；Release 资产建议命名 `LKlight-windows-x64.zip` |

LKlight 源码仓库不直接提交二进制文件。预编译文件应作为 GitHub Release assets 分发，并与 `LICENSE`、`NOTICE`、`README.md` 一同打包，以满足 GPL 源码可得性和署名要求。

### 5.5 编译与运行

```bash
# 编译发行版
cargo build --release

# 运行测试
cargo test --lib    # 29/29 单元测试

# 创建单 swarm 输入
./target/release/LKlight setup tests/1azp/1azp_receptor.pdb tests/1azp/1azp_ligand.pdb -s 1 -g 200

# 运行对接（单 swarm）
./target/release/LKlight run setup.json initial_positions_0.dat 100 pydock

# 完整流水线（所有 swarms 并行）
./target/release/LKlight pipeline receptor.pdb ligand.pdb pydock --threads 8
```

---

## 6. Discussion

### 6.1 主要发现

本工作的核心发现可以概括为两点：

**发现一：Rust 基线版本存在系统性缺陷，无法用于生产环境。** 原始 `lightdock-rust` 二进制文件中，DFIRE/DFIRE2/DDNA 评分函数因外部参数文件缺失而在运行时崩溃，ANM 支持存在 stride 计算错误，多处代码在遇到边界条件（未知残基、ANM 原子数不匹配）时直接 panic 而非优雅降级。这些问题使得原始 Rust 版本实际上不可用于标准 LightDock 工作流。本工作通过全面 Bug 修复、公开单元测试和开发阶段综合数值对比，产出了功能完整的 Rust LightDock 实现。

**发现二：rayon 并行化 + SIMD 友好热路径使 LKlight 全面超越 Python，且覆盖所有主要评分函数。** 通过将受体原子外循环并行化，并把热路径重构为编译器易优化的连续数组与简单内循环，LKlight 全部四个测试场景均超越 Python：pydock **3.0×**、dna+ANM **16.5×**、cpydock **19.2×**、dfire **25.5×**。本机 benchmark 可进一步启用 native CPU 优化观察硬件上限，但公开 Release 二进制保持便携 baseline。

### 6.2 G3/G4 经验教训

空间网格优化是计算结构生物学中的经典加速技术。然而本工作表明，**其有效性高度依赖于截断距离与系统规模的比值**：

- 当 $r_{\text{cut}} / r_{\text{protein}} < 0.3$ 时，空间稀疏，网格有效（如 sd.rs 9Å）
- 当 $r_{\text{cut}} / r_{\text{protein}} > 1.0$ 时（如 pydock 30Å 截断 / 小蛋白配体），大多数原子对均在截断内，查询开销超过计算减少

未来可采用更小的格点（例如 5Å，±6格=13³=2197 — 更差；或 VDW_CUTOFF 10Å 但仅用于 VDW 热路径，ELEC 保持 O(N²)），或改用 SIMD 向量化同时计算 8 个原子对的距离。

### 6.3 局限性

1. **单 swarm 基准**：本基准在单 swarm 维度测量，G2（GSO 运动阶段并行）加速在多 swarm 全流水线场景下更显著，待补充全流水线基准
2. **FP 累加顺序与 Python 的微小差异**：并行化改变了浮点累加顺序，尤其对 dfire/pydock 中的累加和项，可能产生 ~2×10⁻¹³ 级别的数值差异（测试改为 approx 容差、结果内容不变）
3. **更大系统未测试**：当前基准仅针对 1PPE（1615 受体原子）和 1AZP，较大系统（>5000 原子）的 Roofline 行为待验证

---

## 7. Conclusion

本文介绍了 **LKlight v1.0**，一个功能完整、经全面测试验证的 LightDock 分子对接引擎高性能 Rust 实现。主要成果：

1. **修复了原 Rust 基线版本的 4 个系统性 Bug**（DFIRE 崩溃、ANM stride 错误、未知残基 panic、atom_count 断言），产出第一个可用于生产的 Rust LightDock 实现；
2. **公开仓库通过 29/29 单元测试，开发阶段通过 160/160 综合数值对比**，全部 12 类评分函数与 Python 参考实现数值吻合；
3. **全部测试场景全面超越 Python**：pydock **3.0×**、dna+ANM **16.5×**、cpydock **19.2×**、dfire **25.5×**；vs Rust-orig：pydock **26.5×**、dna **307×**、cpydock **163×**（dfire Rust-orig 崩溃，LKlight 为可用 Rust 实现）；
4. **揭示了 G3/G4 HashMap 网格的性能回退根因**（343 次查询开销 >> O(N²) 直接遍历）并通过 H1（rayon 并行化）+ H2（SIMD 友好热路径）完成修复，实现突破性性能提升；
5. **证实 Rust rayon + 编译器友好数据布局可超越 Python NumPy 隐式 SIMD**，无需手写汇编，为 Rust 科学计算实践提供参考。

### 后续工作

| 优先级 | 方向 | 预期收益 |
|--------|------|---------|
| **高** | 全流水线多 swarm 基准（N=400 swarms） | 量化 G2 swarm 并行化实际效果 |
| **中** | 配体坐标向量化变换（rot_mat × coords broadcast） | ANM 坐标更新加速 |
| **中** | 更大测试系统（>5000 原子）基准 | 验证线性/超线性扩展 |
| **低** | GPU 后端（wgpu / CUDA）评分内核 | 超大系统批量筛选 |

---

## 8. GPL-3.0 License Compliance and Attribution

### 8.1 License Status

LKlight is a **derivative work** of LightDock (Python, GPL-3.0-or-later) and its companion Rust implementation `lightdock-rust` (GPL-3.0). Pursuant to GPL-3.0 Section 5, all modifications and extensions distributed as LKlight are released under the **same GPL-3.0-or-later license**.

### 8.2 Is Renaming Compliant?

**Yes.** GPL-3.0 imposes no restriction on the name of a derivative work. The license text (GPL-3.0 §5, §6) requires only that:

| Requirement | LKlight Status |
|-------------|----------------|
| Retain the GPL-3.0 license | `LICENSE` file ✓ |
| Preserve copyright notices | `NOTICE` file ✓ |
| Document significant changes | `CHANGELOG.md` ✓ |
| Make source code available | GitHub publication ✓ |
| Not impose additional restrictions | No additional restrictions ✓ |

The FSF explicitly confirms that GPL permits renaming: *"You may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL."* The name change from `lightdock-rust` to `LKlight` is fully compliant provided all above conditions are met.

### 8.3 What Is NOT Permitted

- Distributing LKlight binaries **without** making the source code available
- Incorporating LKlight into a **proprietary closed-source product** without complying with GPL copyleft obligations
- Removing or obscuring the original LightDock copyright notices
- Claiming that LKlight is the official LightDock release or implying endorsement by the original authors

### 8.4 GitHub Publication Checklist

To publish LKlight on GitHub in full GPL compliance, the following files are required:

| File | Purpose | Status |
|------|---------|--------|
| `LICENSE` | Full GPL-3.0 license text | ✓ Present |
| `NOTICE` | Copyright attribution to original LightDock/lightdock-rust authors | ✓ Present |
| `CHANGELOG.md` | Record of significant changes from upstream | ✓ Present |
| `CONTRIBUTING.md` | Contribution guidelines | ✓ Present |
| `README.md` | Project description with license badge | ✓ Present |
| `Cargo.toml` | Correct `license = "GPL-3.0-or-later"` field | ✓ Present |
| Source code | All Rust source files | ✓ `src/` directory |
| `Cargo.lock` | Reproducible binary build | ✓ Present |
| `.github/workflows/rust.yml` | Cross-platform CI | ✓ Present |
| `RELEASE.md` | Release and binary asset guidance | ✓ Present |
| `.gitattributes` | Preserve binary parameter files | ✓ Present |

> **Note:** The `Cargo.lock` file should be committed for binary executables (as recommended by Cargo). The `target/` build directory should remain in `.gitignore`.

### 8.5 Current Publication Status

The public source repository is:

```text
https://github.com/LK-Studio1128/LKlight
```

Generated build artifacts (`target/`, `dist/`, `.DS_Store`, backup files) are intentionally excluded from Git. Pre-built macOS/Linux/Windows binaries should be uploaded as GitHub Release assets rather than committed to the source tree.

---

## Acknowledgements

LKlight builds upon the intellectual and engineering foundations of LightDock, developed by Brian Jiménez-García, Jorge Roel-Touris, and collaborators at the Life Sciences Department, Barcelona Supercomputing Center (BSC), Spain. We gratefully acknowledge their development of the original LightDock framework [1,2] and the `lightdock-rust` Rust baseline, both made freely available under GPL-3.0. The `rayon` parallel data processing library [13] and the Rust compiler's LLVM backend are essential enabling infrastructure for the performance results reported here. Benchmark PDB structures 1PPE (trypsin-BPTI complex) and 1AZP (protein-DNA complex) were retrieved from the RCSB Protein Data Bank [14].

---

## References

1. Jiménez-García B, Roel-Touris J, Romero-Durana M, Vidal M, Jiménez-González D, Fernández-Recio J. **LightDock: a new multi-scale approach to protein–protein docking.** *Bioinformatics.* 2018;34(1):49–55. doi:10.1093/bioinformatics/btx555

2. Roel-Touris J, Bonvin AMJJ, Jiménez-García B. **LightDock goes information-driven.** *Bioinformatics.* 2020;36(3):950–952. doi:10.1093/bioinformatics/btz642

3. Krishnanand KN, Ghose D. **Glowworm swarm optimization for simultaneous capture of multiple local optima of multimodal functions.** *Swarm Intelligence.* 2009;3(2):87–124. doi:10.1007/s11721-008-0021-5

4. Zhou H, Zhou Y. **Distance-scaled, finite ideal-gas reference state improves structure-derived potentials of mean force for structure selection and stability prediction.** *Protein Sci.* 2002;11(11):2714–2726. doi:10.1110/ps.0217002

5. Yang Y, Zhou Y. **Specific interactions for ab initio folding of protein terminal regions with secondary structures.** *Proteins.* 2008;72(2):793–803. doi:10.1002/prot.21968

6. Cheng TM, Blundell TL, Fernandez-Recio J. **pyDock: electrostatics and desolvation for effective scoring of rigid-body protein–protein docking.** *Proteins.* 2007;68(2):503–515. doi:10.1002/prot.21419

7. Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin O, Bahar I. **Anisotropy of fluctuation dynamics of proteins with an elastic network model.** *Biophys J.* 2001;80(1):505–515. doi:10.1016/S0006-3495(01)76033-X

8. Mintseris J, Pierce B, Wiehe K, Anderson R, Chen R, Weng Z. **Integrating statistical pair potentials into protein complex prediction.** *Proteins.* 2007;69(3):511–520. doi:10.1002/prot.21502

9. Feliu E, Oliva B. **How different from random are docking predictions when scoring is poor?** *J Chem Inf Model.* 2010;50(12):2153–2160. doi:10.1021/ci100369y

10. Miyazawa S, Jernigan RL. **Residue-residue potentials with a favorable contact pair term and an unfavorable high packing density term, for simulation and threading.** *J Mol Biol.* 1996;256(3):623–644. doi:10.1006/jmbi.1996.0114

11. Ravikant DVS, Elber R. **PIE—efficient filters and coarse grained potentials for unbound protein-protein docking.** *Proteins.* 2010;78(2):400–419. doi:10.1002/prot.22566

12. Matsakis ND, Klock FS II. **The Rust language.** *ACM SIGAda Ada Letters.* 2014;34(3):103–104. doi:10.1145/2692956.2663188

13. Stone J, et al. **Rayon: A data parallelism library for Rust.** https://github.com/rayon-rs/rayon (accessed 2025).

14. Berman HM, Westbrook J, Feng Z, et al. **The Protein Data Bank.** *Nucleic Acids Res.* 2000;28(1):235–242. doi:10.1093/nar/28.1.235

15. Van Zundert GCP, Rodrigues JPGLM, Trellet M, et al. **The HADDOCK2.2 Web Server: User-Friendly Integrative Modeling of Biomolecular Complexes.** *J Mol Biol.* 2016;428(4):720–725. doi:10.1016/j.jmb.2015.09.014

16. Trott O, Olson AJ. **AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading.** *J Comput Chem.* 2010;31(2):455–461. doi:10.1002/jcc.21334

17. Ester M, Kriegel H-P, Sander J, Xu X. **A density-based algorithm for discovering clusters in large spatial databases with noise.** *Proc KDD.* 1996:226–231. (DBSCAN, basis for LightDock cluster subcommand)

---

## Appendix A — Optimization Timeline

| Session | 优化项 | 状态 |
|---------|--------|------|
| 3 | glowworm.rs Vec→[f64;3]；DFIRE/DFIRE2/DDNA thread-local HashMap 复用 | ✓ |
| 3 | pisa.rs + tobi.rs 空间索引（O(N²)→O(N)）；tobi.rs 消除 sqrt | ✓ |
| 4 | Bug 修复（Fix 1-4）；160/160 综合测试通过 | ✓ |
| 5 | F1 sqrt_vdw_charges；F2 sd.rs 9Å 网格；F3 BufWriter；F4 qt [f64;3] | ✓ |
| 6 | G1 pos/rot scratch；G2 movement 并行；G3/G4 HashMap 网格（产生回退） | G3/G4 已撤回 |
| 7 | **H1** pydock/dna/cpydock rayon 并行外循环；**H2** SIMD 友好热路径与可选 native benchmark 构建 | ✓ **pydock 3.0×Py, 26.5×Orig** |
| 8 | **I1** dfire rayon并行+移除 HashMap；**I2** dfire2 同；**I3** sd.rs 并行化 | ✓ **dfire 25.5×Py** (935ms→33ms) |

## Appendix B — Test Command Reference

```bash
# 单元测试
cargo test --lib 2>&1 | tail -5

# 性能基准
bash benchmark.sh

# 轻量夹具 smoke test（1AZP pydock，swarm 0）
mkdir -p demo-1azp
cd demo-1azp
../target/release/LKlight setup ../tests/1azp/1azp_receptor.pdb \
    ../tests/1azp/1azp_ligand.pdb -s 1 -g 200
../target/release/LKlight run setup.json initial_positions_0.dat 100 pydock
```
