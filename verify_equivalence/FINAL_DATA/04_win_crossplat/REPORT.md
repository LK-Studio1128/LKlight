# LKlight Windows x86-64 跨平台验证报告

**日期**: 2026-08-28
**验证人**: LKlight 验证管道（自动）

## 环境对比

| 项 | macOS 端 | Windows 端 |
|---|---|---|
| 主机 | Apple Silicon (darwin) | 117.50.173.234 (Win Server 2022, build 20348) |
| CPU | Apple Silicon arm64 | Intel Xeon Gold 5118 @ 2.30GHz (AMD64) |
| 二进制 | LKlight (arm64, 2026-08-27 23:48) | LKlight.exe (交叉编译 x86_64-pc-windows-gnu, PE32+, 2026-08-28 03:10) |
| 编译器 | rustc 1.94.1 | rustc 1.94.1 + x86_64-w64-mingw32-gcc |

## 验证设计

- **体系**: BM5 2X9A（蛋白-蛋白复合物，fastdfire 评分）
- **参数**: swarms=10 × glowworms=50 × 100 steps，seed=42
- **输入**: 官方 lightdock3_setup.py 表面点初始位置（macOS 端生成后上传，两平台共享同一初始位置与 setup.json）
- **对照**: 同一输入在 macOS 端 30 例批量验证中的 2X9A 运行结果

## 结果

**110/110 检查点文件字节级完全一致**（10 swarm × 11 个检查点 [步 1,10,20,...,100] × 每文件 50 条 glowworm 状态记录，每条含 7 维坐标 + luciferin 等 11 个数值字段，9 位小数精度）。

- 数值最大偏差：**0.0（逐位一致）**
- glowworm 数量/顺序：完全一致
- 跨平台确定性成立：GSO 轨迹在 macOS arm64 与 Windows x86-64 上逐位相同

## 结论

Windows x86-64 版 LKlight.exe 与 macOS arm64 版在相同输入下产生**逐位相同**的对接轨迹，证实：
1. 修复后的代码（含 4 项官方基线缺陷修复 + 6 项移植偏差修复）在 Windows 平台行为一致；
2. 精度等价结论（30 例 BM5 验证，平均 top-N 差 ≤0.02）可推广到 Windows 平台；
3. 论文"跨平台单二进制"声明有了直接的数值证据支撑。

## 文件清单

- swarm_0..9/ — Windows 端 GSO 输出（11 个检查点文件/swarm）
- score_timings.csv, score_timings_ps.csv — 此前性能计时（2 核 Xeon 限制下仅供参考）
- numeric_win.csv — 此前数值验证记录
