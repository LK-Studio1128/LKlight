# LKlight v1.1.0

[![CI](https://github.com/LK-Studio1128/LKlight/actions/workflows/rust.yml/badge.svg)](https://github.com/LK-Studio1128/LKlight/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-stable-orange.svg)](https://rustup.rs)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22150512.svg)](https://doi.org/10.5281/zenodo.22150512)

**LKlight** 是 [LightDock](https://lightdock.org) 分子对接引擎的高性能 Rust 再实现（GPL-3.0 衍生作品）。本项目基于 LightDock 的 GSO（Glowworm Swarm Optimization）分子对接思想和上游 [`lightdock-rust`](https://github.com/lightdock/lightdock-rust) Rust 基线继续开发，目标是在保留 LightDock 方法体系与可复现性的前提下，提供更快、更稳定、更易分发的单文件命令行引擎。

> 基于 LightDock GSO 算法 | 12 类评分函数 | rayon 并行 | SIMD 友好 | 全参数内嵌、单文件分发

---

## 项目来源与 LKlight 的工作

LKlight 不是从零开始的新算法，而是对 LightDock / lightdock-rust 的工程化增强和性能优化版本。我们保留并明确标注原始 LightDock 作者、论文和 GPL 许可，同时将上游 Rust 原型扩展为可直接用于自动化流程和跨平台发布的高性能 docking engine。

### 源项目

| 项目 | 说明 |
|------|------|
| [LightDock](https://lightdock.org) | 原始 Python 分子对接框架，提出基于 GSO 的多尺度蛋白质-蛋白质 / 蛋白质-DNA 对接流程 |
| [`lightdock-rust`](https://github.com/lightdock/lightdock-rust) | LightDock 官方 Rust 基线实现，提供核心运行框架和部分评分函数 |
| **LKlight** | 在上游 Rust 基线基础上进行修复、扩展、并行化、数据内嵌和跨平台二进制发布 |

### 我们完成的主要工作

| 类别 | 工作内容 |
|------|----------|
| 评分函数扩展 | 支持 12 类评分函数、13 个命令行方法名：`dfire`、`fastdfire`、`dfire2`、`dna`、`mj3h`、`pydock`、`cpydock`、`sd`、`vdw`、`pisa`、`sipper`、`tobi`、`ddna` |
| 稳定性修复 | 修复 DFIRE/DFIRE2/dDNA 外部参数缺失导致的运行时崩溃；修复 ANM stride、非标准残基、ANM atom-count assertion 等上游基线问题 |
| 数据内嵌 | 将 DFIRE、DFIRE2、dDNA、MJ、PISA 等关键参数随源码/二进制分发，减少运行时外部文件依赖 |
| 性能优化 | 使用 rayon 并行化受体原子外循环，重构热路径为 SIMD 友好形式，复用 thread-local scratch buffer，减少堆分配 |
| 算法工程化 | 为 SD/PISA/TOBI 等短截断或接触势路径加入空间索引/剪枝策略，降低无效 pair 计算 |
| CLI 整合 | 提供统一 `LKlight` 单二进制入口，覆盖 setup、run、rank、cluster、generate、score、pipeline、trajectory 等常用流程 |
| 跨平台发布 | 提供 macOS、Linux、Windows 构建脚本；预编译二进制可作为 GitHub Release 资产独立分发 |
| 开源合规 | 保留 GPL-3.0-or-later、NOTICE、CHANGELOG、CONTRIBUTING，并明确上游来源和论文引用 |

### 适用场景

- 蛋白质-蛋白质 docking
- 蛋白质-DNA docking
- 含 ANM 柔性模式的 docking
- 需要批处理、自动化脚本或单文件二进制部署的本地 docking 工作流
- 需要对 LightDock Python 版本或 lightdock-rust 基线进行性能对比的研究/工程场景

---

## 性能基准

测试平台：macOS arm64（Apple Silicon），swarm_0，200 glowworms，100 步，3 次均值

| 场景 | Python | Rust-orig | **LKlight** | LKlight/Py | LKlight/Orig |
|------|--------|-----------|-------------|-----------|-------------|
| 1PPE pydock | 858 ms | 7,693 ms | **290 ms** | **3.0×** | **26.5×** |
| 1PPE dfire | 840 ms | CRASH ¹ | **33 ms** | **25.5×** | — |
| 1AZP dna+ANM | 760 ms | 14,142 ms | **46 ms** | **16.5×** | **307×** |
| 1PPE cpydock | 844 ms | 7,158 ms | **44 ms** | **19.2×** | **163×** |

> ¹ Rust-orig dfire 因外部参数文件缺失在运行时崩溃；LKlight 将参数完整内嵌到二进制中。

---

## 支持的评分函数

| 方法名 | 描述 |
|--------|------|
| `dfire` | DFIRE 统计势（参数内嵌） |
| `fastdfire` | `dfire` 的兼容别名 |
| `dfire2` | DFIRE2 统计势 |
| `dna` | 蛋白质-DNA 评分 |
| `pydock` | PyDock 静电+VdW |
| `cpydock` | cpyDOCK（含去溶剂化） |
| `sd` | CHARMM 力场（含切换函数） |
| `vdw` | 纯 VdW |
| `mj3h` | MJ 接触势 |
| `pisa` | PISA 溶剂化势 |
| `sipper` | SIPPER 接触矩阵 |
| `tobi` | TOBI 接触势 |
| `ddna` | DNA 专用势 |

---

## 安装

需要 [Rust 工具链](https://rustup.rs/)（推荐 stable 1.75+）。

```bash
git clone https://github.com/LK-Studio1128/LKlight.git
cd LKlight
cargo build --release
# 产物：target/release/LKlight
```

macOS 一键打包（输出 `dist/LKlight-mac/LKlight`）：

```bash
bash build_mac.sh
```

Linux 静态二进制（musl，适用所有发行版）：

```bash
bash build_linux.sh
```

---

## 快速开始

### 设置对接系统

```bash
LKlight setup receptor.pdb ligand.pdb -s 25 -g 200
# -s 25   : 25 个 swarm
# -g 200  : 每 swarm 200 个 glowworm
# --anm   : 开启 ANM 柔性（可选）
```

### 运行优化

```bash
LKlight run setup.json initial_positions_0.dat 100 pydock
#                                               ^^^  步数
```

### 完整子命令

```
LKlight setup      <rec.pdb> <lig.pdb> [-s N] [-g N] [--anm] [--restraints F]
LKlight run        <setup.json> <initial_positions.dat> <steps> <method>
LKlight rank       <num_swarms> <steps> [--filter-clusters]
LKlight rank_swarm <num_swarms> <steps>
LKlight cluster    <gso_output.out> [--cutoff 4.0]
LKlight top        <ranking_file> <N>
LKlight filter     <ranking_file> <restraints.list>
LKlight generate   <rec.pdb> <lig.pdb> <gso.out> <N>
LKlight score      <rec.pdb> <lig.pdb> <method> [--tx X --ty Y --tz Z]
LKlight diameter   <pdb_file>
LKlight gso_to_csv <ranking_file> <output.csv>
LKlight move_anm   <pdb> <n_modes> <n_confs>
LKlight map_contacts <rec.pdb> <lig.pdb> <gso_file>
LKlight trajectory <rec.pdb> <lig.pdb> <swarm_id> <glowworm_id> <steps>
LKlight reference_points <pdb_file> [--save]
LKlight pipeline   <rec.pdb> <lig.pdb> <method> [--threads N]
```

---

## 使用示例

### 1AZP — 蛋白质-配体/蛋白质-DNA 测试夹具（pydock）

```bash
mkdir -p demo-1azp
cd demo-1azp
../target/release/LKlight setup ../tests/1azp/1azp_receptor.pdb ../tests/1azp/1azp_ligand.pdb -s 1 -g 200
../target/release/LKlight run setup.json initial_positions_0.dat 100 pydock
```

```
Swarm ID 0 | pydock | 200 glowworms | 100 steps
Writing to swarm dir "swarm_0"
Done.
```

### 1AZP — ANM 柔性示例

```bash
mkdir -p demo-1azp-anm
cd demo-1azp-anm
../target/release/LKlight setup ../tests/1azp/1azp_receptor.pdb ../tests/1azp/1azp_ligand.pdb -s 1 -g 200 --anm
../target/release/LKlight run setup.json initial_positions_0.dat 100 dna
```

```
Swarm ID 0 | dna+ANM | 200 glowworms | 100 steps
Writing to swarm dir "swarm_0"
Done.
```

### DFIRE 统计势评分

```bash
../target/release/LKlight score ../tests/1azp/1azp_receptor.pdb ../tests/1azp/1azp_ligand.pdb dfire
```

```
Score: <value>
```

---

## 测试

```bash
cargo test --lib          # 29 单元测试
bash benchmark.sh         # 三版本性能对比（需 Python lightdock3）
```

---

## 主要优化

| 版本 | 优化项 |
|------|--------|
| v1.0 | rayon 并行受体原子外循环（pydock/dna/cpydock/dfire/dfire2/sd） |
| v1.0 | SIMD 友好的连续数组与内层循环设计；发布配置使用便携 CPU baseline，benchmark 可按需启用 native 优化 |
| v1.0 | DFIRE/DFIRE2 参数内嵌二进制（消除外部文件依赖） |
| v1.0 | thread-local 坐标/界面缓冲区复用（消除每步堆分配） |
| v1.0 | sd.rs 9Å 空间网格剪枝（真正稀疏短截断） |
| v1.0 | pisa/tobi 空间索引（O(N²)→O(N)） |

---

## 贡献

欢迎提交 Issue 和 Pull Request！详见 [CONTRIBUTING.md](CONTRIBUTING.md)。

---

## License & Attribution

LKlight 是 [LightDock](https://lightdock.org) 的 **GPL-3.0 衍生作品**。

- 许可证：**GPL-3.0-or-later**（继承自 LightDock）
- 完整许可证文本见 [LICENSE](LICENSE)
- 原始版权归属见 [NOTICE](NOTICE)
- 修改记录见 [CHANGELOG.md](CHANGELOG.md)

**必须保留出处：** Brian Jiménez-García 等人在巴塞罗那超级计算中心开发了原始 LightDock 框架。
请引用： Jiménez-García B, et al. *Bioinformatics* 2018;34(1):49–55. doi:[10.1093/bioinformatics/btx555](https://doi.org/10.1093/bioinformatics/btx555)

