# HPC 升档执行包（hpc_bundle）

> 本包是 `方案二` 两条"升档硬条件"的**自包含、可移植执行单元**，由计算沙箱制备、
> 供在 **HPC（SCNet/超算，开放网络 + vina 二进制）** 上真实执行。
> 沙箱本身受网络/编译约束无法跑分，故**真实结果只能在此包于 HPC 执行后产生**。

## 包含内容（已镜像仓库目录结构）

| 路径 | 作用 |
|---|---|
| `00_pipeline/config.py` | 中央配置（ROOT/RESULTS/THRESH 等；相对路径由此解析） |
| `01_data/lincs/gene_info.txt.gz` | LINCS L1000 landmark 基因映射（step0b 构建签名所需） |
| `01_data/lincs/control_sig_ids.txt` | 7 个阳性对照药在 L1000 的 **112** 个 sig_id |
| `02_scripts/step0_lincs.py` | compute_cs / build_query / load_gene_symbol_map 等（step0b 复用） |
| `02_scripts/step0b_lincs_g1_controls.py` | Step 0b G1 重算主程序 |
| `02_scripts/download_lincs_lvl3_hpc.sh` | HPC 下载完整 LINCS Level3（~5GB，断点续传） |
| `03_results/step1_activation_up_full.csv` `..._down_full.csv` | NP 激活签名（step0b 构建 landmark 索引所需） |
| `06_docking/run_vina.sh` | 30 条 vina 调用的批量脚本 |
| `06_docking/run_vina.slurm` | SLURM 调度头（8 核 / 12h / 16G） |
| `06_docking/manifest.csv` | 30 任务清单 |
| `06_docking/rec_pdbqt/*.pdbqt`（6） | 受体 PDBQT（已剔除非标准残基） |
| `06_docking/lig_pdbqt/*.pdbqt`（5） | 配体 PDBQT（含完整扭转树） |

## 执行

```bash
tar -xzf hpc_bundle.tar.gz && cd hpc_bundle
module load vina            # 或 conda activate 含 vina 的环境
bash RUNME.sh               # 自动：(a) 下载+重算 G1  (b) 跑 Vina 对接
```

或分步：
```bash
# (a) 仅 G1
bash 02_scripts/download_lincs_lvl3_hpc.sh 01_data/lincs
python3 02_scripts/step0b_lincs_g1_controls.py \
        --gctx 01_data/lincs/GSE92742_LVL3.gctx \
        --controls 01_data/lincs/control_sig_ids.txt

# (b) 仅对接
bash 06_docking/run_vina.sh          # 或 sbatch 06_docking/run_vina.slurm
```

## 产物与回填

- **G1**：`03_results/step0b_summary.json` → `g1_pass` / `mwu_one_sided_p` / 逐药 CS。
  - 回填规则：仅当 `g1_pass=true` 时，把 Step 7/8 重定位由"探索性"升为"闸门验证通过"，并改写 v4 §3.7'/闸门表。
- **Vina**：`06_docking/results/*.log` → 每任务 9 个 pose 的 binding affinity（kcal/mol）。
  - 回填规则：仅引用真实 `*.log` 中的数值；若最佳 pose ≤ −7 kcal/mol 且优于随机/decoy，则把对接线由"已备料"升为"硬阳性"，并触发向发现型 Q1 的升档重评（见 v4 §6 期刊纪律）。

## 资源提示

- **Step 0b 内存**：完整 Level3 矩阵约 1.3M 签名 × 978（float32 ≈ 5GB / float64 ≈ 10GB）。
  脚本以 `mat = f[...][:]` 一次性读入内存；请在 **≥16GB 内存**节点运行（SLURM 已申请 16G）。
  如内存紧张，可改为分块读取（修改 `step0b_lincs_g1_controls.py` 的 `run()`），但默认实现需足够 RAM。
- **Vina 时长**：30 任务 × `--exhaustiveness 32` 在 8 核约数分钟~数十分钟（视结构大小）；SLURM 申请 12h 充裕。

## 诚实边界（不可突破）

在 HPC 真实产物回传前，**不得**在手稿中引用任何 G1 判定或对接 kcal/mol 数字。
本包只负责"把已备料的计算真实地跑出来"；任何升档结论都须由这里产出的数字支撑，严禁编造。
