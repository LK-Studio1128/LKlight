# LKlight 全面测试与 Mac GPU(Metal) 集成报告 — 2026-09-06

测试主机: Mac mini M4 (10-core, 32 GB, Metal 4) · 测试负载: 1AZP, dna, seed 324324,
1 swarm × 1000 glowworms × 1000 steps (与 5090 验收 / 论文表 4 完全同参数)。

## 1. 同架构便携性(换机即插即用) — 通过

| 二进制 | 动态依赖 (otool -L) | 签名 | 结论 |
|---|---|---|---|
| exact mac-arm64 v1.2.0 | 仅 /usr/lib/libSystem, libiconv | arm64 adhoc linker-signed | 任意 Apple Silicon 拷贝即用 |
| grid mac-arm64 v1.2.0 | 仅 /usr/lib/libSystem, libiconv | arm64 adhoc linker-signed | 同上 |
| **LKlight-metal-arm64 (新)** | /System Metal+Foundation+CoreFoundation, /usr/lib objc/libSystem/libiconv | arm64 adhoc linker-signed | 同上(Metal 为 macOS 系统框架,所有 Apple Silicon 具备) |

- 无任何 homebrew / @rpath / 绝对路径 / 机器绑定依赖; 引擎只按相对 cwd 读
  setup.json 与 PDB, 换目录/换机无需重装。
- 佐证: 同一引擎二进制在三个不同位置(原构建目录、洁净拷贝目录、以及跨机器
  Windows 5090 / Linux 3080Ti / Mac M4)全量跑出的 best 能量逐位一致(§2)。
- 附注: x86_64 目标在 .cargo/config.toml 中固定 target-cpu=x86-64 保证 x86 换机兼容;
  aarch64 使用默认 (Apple Silicon 基线), 同样与构建机无关。

## 2. 全量回归矩阵(洁净目录复跑, 空闲机独占)

| 引擎 | 位置 | 墙钟 | Best 能量 | 与基线比对 |
|---|---|---|---|---|
| exact | ht1000_mac/exact (02:19) | 247.83 s | -7251.9026691900 | = 5090 逐位 |
| grid | ht1000_mac/grid (02:24) | 57.56 s | -7254.9522016100 | = 5090 逐位 |
| exact(复跑,洁净拷贝) | ht1000_pristine/exact | 219.29 s | -7251.9026691900 | 逐位一致 ✅ |
| grid(复跑,洁净拷贝) | ht1000_pristine/grid | 47.10 s | -7254.9522016100 | 逐位一致 ✅ |
| gpu-repo CPU 回退 | ht1000_gpurepo_cpu | 41.34 s | -7254.9522016100 | = 独立 grid 引擎逐位 ✅ |
| **metal (新,集成版)** | ht1000_metal | **6.22 s** | **-7254.9516108000** | = 5090 CUDA 引擎 -7254.9516108 逐位 ✅ |

交叉验证要点:
- 跨平台跨位置逐位一致: Windows RTX5090 / Linux RTX3080Ti / macOS M4 三平台
  exact & grid best 完全相同 → 引擎确定性与代码同源性强。
- gpu 仓库默认 CPU 构建(即 metal 不可用时的回退路径)与独立 grid 引擎 best
  逐位一致 → 任何加速器失效时静默回退 CPU, 结果可预期。

## 3. Mac GPU (Metal) 引擎 — 正式接入完成

新增(均在 LKlight-GPU 仓库, `--features metal` 门控, 默认构建不受影响):
- `src/metal/lk_metal.m` — ObjC shim + 内嵌 MSL kernel (POC v4 终版逐字移植:
  per-pose threadgroup ×128, Kahan 补偿累加+two-sum 树归约, 远场 φ 三线性 + 近场
  27-cell cell-list)。静态库经 build.rs 由 clang 编译链接 (Metal/Foundation)。
- `src/metal_score.rs` — Rust 绑定, 镜像 `gpu_score.rs`: OnceLock 持久 ctx、
  `metal_available()`、`batch_energy_metal_scores()`。CPU f64 刚体变换留在 host
  (绕开 Apple GPU 无 fp64), 共享 MTLBuffer 零拷贝上传。
- `src/dna.rs` — `lig_metal` OnceLock ctx + `supports_batch`/`batch_energy`
  cfg 分派: metal → cuda → CPU grid 三级优雅回退。

引擎级实测 (1000×1000 全量, 含每步文件 I/O 与群算法):

| 指标 | 数值 |
|---|---|
| Metal 引擎总墙钟 | **6.22 s** (6.2 ms/步 × 1000) |
| vs 同场 CPU grid (47.10 s) | **7.6×** |
| vs 同场 CPU exact (219.29 s) | **35.3×** |
| best 能量 | -7254.95161080 (= 5090 CUDA 终值, 逐位) |
| 数值档位 | f32 量化 ~1e-4 相对 (与 CUDA 版同档, 论文脚注口径一致) |

## 4. 结论

- 同架构换机使用: 满足。三引擎二进制纯系统依赖 + adhoc 签名 + 位置无关。
- Mac GPU 加速: 满足。Metal 引擎跑通 1000×1000 验收, 端到端 7.6× (grid)/ 35.3×
  (exact), best 与 CUDA 版逐位一致; 无 GPU 时自动回退 CPU 且结果与 grid 引擎逐位一致。
- 待用户拍板: 版本发布 (v1.2.1: metal feature 四平台 + mac metal 资产) 与论文
  §3.2/表 4 补 Mac GPU 行。

复现命令:
- metal 引擎: `cargo build --release --features metal` → 运行 `LKlight run setup.json initial_positions_0.dat 1000 dna`
- CPU: 去掉 `--features metal` 即可 (回退路径与 grid 引擎同结果)。
