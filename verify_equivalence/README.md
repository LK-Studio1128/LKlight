# LKlight / LightDock Accuracy-Equivalence Benchmark

验证 **LKlight** 与 **Python LightDock** 在标准对接任务上的精度等价性 —— 这是投稿
《LKlight: filling the functional and evidential gaps of the official LightDock Rust
engine》所必需的 P0 实验。

## 原理

对 Protein-Protein Docking Benchmark 5（`lightdock_bm5`）中每个复合物：

1. **共享输入**：用 `LKlight setup` 生成 `setup.json` + `initial_positions_*.dat`
   （两版引擎共用同一套初始位置，保证公平）。
2. **LKlight**：`LKlight run`（每 swarm）→ `LKlight generate`（生成全部构象 PDB）。
3. **Python LightDock**（二选一）：
   - `--py-mode run`：本地运行 `lightdock3` + `lgd_generate_conformations.py`；
   - `--py-mode official`：直接解析官方 `results/*.list`（官方已跑好的 Python 版
     iRMSD/LRMSD/Fnat/score），零算力成本，与官方 2020 论文数据同源。
4. **评估**：按评分降序取 top-N 构象，与 native 结构比对，按 CAPRI 标准分级
   （Acceptable/Medium/High），成功率 = top-N 内含 ≥1 个可接受模型的比例。

论文中报告的成功率口径与 LightDock 2018/2020 官方论文一致（success rate @ top-N）。

## 环境准备

```bash
# 1) 构建 LKlight
cd LKlight && cargo build --release

# 2) 获取官方 benchmark 数据（~小仓库，浅克隆即可）
git clone --depth 1 https://github.com/lightdock/lightdock_bm5.git

# 3a) 需要自己跑 Python 版时（--py-mode run）：
pip install lightdock           # 提供 lightdock3 / lgd_generate_conformations.py
# 3b) 只想用官方基线时（--py-mode official）：无需安装 Python 版

# 4) 必需：numpy（L-RMSD 的配体 Kabsch 对齐依赖它；缺失时评估精度不可靠）
#    本机 macOS 可直接用自带 numpy 的解释器，如 /opt/homebrew/bin/python3
pip install numpy matplotlib     # matplotlib 仅用于出图（可选）
```

## 快速开始

```bash
cd verify_equivalence

# 推荐（零算力，用官方 Python 结果做基线；LKlight 跑 20 例）
# 注意：用带 numpy 的解释器（如 /opt/homebrew/bin/python3）
python3 run_equivalence.py \
  --bm5 ../bm5_data \
  --lklight ../target/release/LKlight \
  --py-mode official \
  --swarms 10 --glowworms 100 --steps 100 \
  --n-cases 20 --jobs 4

# 完整等价验证（本地同时跑 Python 与 LKlight，同参）
python3 run_equivalence.py \
  --bm5 ../bm5_data \
  --lklight ../target/release/LKlight \
  --py-mode run \
  --py-lightdock $(which lightdock3) \
  --py-generate $(which lgd_generate_conformations.py) \
  --swarms 10 --glowworms 100 --steps 100 --scoring fastdfire \
  --n-cases 20 --jobs 4
```

## 参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `--bm5` | - | `lightdock_bm5` 仓库路径（含 `data/`） |
| `--fetch-bm5 DIR` | - | 自动浅克隆数据仓库到 DIR |
| `--lklight` | `target/release/LKlight` | LKlight 二进制 |
| `--py-mode` | `run` | `run` 本地跑 Python；`official` 用官方结果；`skip` 只跑 LKlight |
| `--swarms` | 10 | 每例 swarm 数（官方 400；同参对比即可，论文注明绝对成功率低于官方高参设置） |
| `--glowworms` | 100 | 每 swarm 萤火虫数（官方 200） |
| `--steps` | 100 | 优化步数（官方 100） |
| `--scoring` | `fastdfire` | 评分函数（官方 BM5 用 fastdfire） |
| `--n-cases` | 20 | 取 BM5 前 N 个复合物（论文建议 20–40） |
| `--systems` | 空 | 显式指定复合物列表，如 `--systems 1PPE 1AZP` |
| `--jobs` | 2 | 复合物级并行数 |

## 输出

- `equivalence_results/summary.tsv`：每例、每引擎的 top-1/5/10/20/50/100 成功率
- `equivalence_results/<case>/per_model_{lk,py}.tsv`：每个 top 模型的
  Fnat / L-RMSD / 受体 RMSD / CAPRI 分级
- `equivalence_results/success_curves.png`：成功率曲线对比图

## 论文里的使用口径（务必阅读）

1. **成功率口径**：与官方论文一致（top-N 内含 ≥1 acceptable 模型的比例）。
   等价性看的是 **LKlight 曲线与 Python 曲线重合**，不是绝对成功率高低。
2. **同参公平**：两版引擎永远使用相同的 swarm/glowworm/step/scoring 和同一套
   `initial_positions`；`--py-mode official` 时官方结果是 400 swarm 高参设置，
   只作参考基线，不可与 LKlight 低参结果直接比绝对值（论文中需注明）。
3. **随机种子**：LKlight 的 `--seed` 保证可复现；Python 版随机种子由 `setup.json`
   的 `seed` 字段控制（`LKlight setup --seed` 会写入）。
4. 若希望完整复现官方成功率（top-10 ≈ 官方水平），可将 `--swarms 400` 升级运行
   （算力：20 例 × 400 swarm × 100 步，多核下需数小时到一天）。

## 常见问题

- **`lgd_generate_conformations.py` 找不到**：lightdock 安装后脚本在 Python bin 目录，
  可用 `python3 -m pip show lightdock` 定位，或直接用 `--py-mode official`。
- **native 结构匹配失败（rec_rmsd=inf）**：多为 `*_A-noh.pdb` / `*_B-noh.pdb`
  命名不符，用 `--systems` 指定具体复合物并检查 `data/<case>/` 下的 PDB 命名。
- **生成构象磁盘占用大**：`glowworms × swarms` 个 PDB/例；完成后可删除
  `equivalence_results/*/lk/swarm_*/lightdock_*.pdb`（`summary.tsv` 已保留结果）。
