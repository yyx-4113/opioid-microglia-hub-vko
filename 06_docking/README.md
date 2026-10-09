# 06_docking — 真实 AutoDock Vina 对接输入包与 HPC 执行说明

> **诚实状态边界（务必先读）**
> 本目录在**计算沙箱（WorkBuddy sandbox）**中完成了**全部输入制备与脚本编写**，并**已真实执行对接与 decoy 特异性对照**：
> 1. `pip install vina` 编译失败（缺 Boost 构建依赖，且本机 Python 3.13 无预编译 wheel）；沙箱代理最初拦截 Vina 二进制。
> 2. **后续取得预编译 AutoDock Vina 1.2.7 Windows 二进制**，得以在本地（沙箱机器）以与 null 同版本重跑 30 个 real-box 任务，并补全 **360 随机盒 + 240 decoy 配体**的特异性对照（见 `decoy_control/`）。
> 3. HPC（SCNet/超算）最初承担真实打分执行（Vina 1.2.3，headline −8.50）；最终 real-box 以本地 Vina 1.2.7 重跑（版本匹配 null），headline = JAK3–amlexanox **−8.53 kcal/mol**，30/30 负值，10/30 ≤ −7。
>
> 因此：**沙箱 = 备料 + 代码/输入校验 + 本地 1.2.7 真实执行（含 decoy 对照）**。decoy/随机盒特异性对照**已完成**（randbox null 精确单尾 p=0.019、decoy null p=0.020，特异性显著但温和）。手稿 §3.10、§5、§6 已据真实产物回填；所有 kcal/mol 均来自真实对接输出（`decoy_control/results_*/`），无编造。

---

## 1. 目录清单（沙箱已制备，已校验非空的产物）

| 类别 | 文件 | 说明 | 完整性校验 |
|---|---|---|---|
| 受体源 PDB | `receptors/*.pdb`（6） | JAK3_6RU9 / IKBKB_4KIK / MALT1_5CQY / TBK1_4DLE / IRAK1_4U2U / MYD88_5NFO（均来自 RCSB，沙箱可直连下载） | 305–1750 KB |
| 受体 PDBQT | `rec_pdbqt/*.pdbqt`（6） | openbabel 3.2.1 + meeko 预处理；已剔除非标准残基（含 JAK3 的 20 个 SEC 硒代半胱） | 158 KB–1.02 MB（**全部非零**） |
| 配体 PDBQT | `lig_pdbqt/*.pdbqt`（5） | tofacitinib / ruxolitinib / phenazone / amlexanox / vorinostat（rdkit 3D 构象 + meeko 写 PDBQT，含 ROOT/H PARENT 扭转树） | 1.3–2.4 KB（**全部非零**） |
| 任务清单 | `manifest.csv` | 6 受体 × 5 配体 = 30 个对接任务；含 `exploratory_box` 标记（MALT1/IRAK1 因无共晶配体/已知口袋，盒子偏大、判为探索性） | 30 行 |
| 批量脚本 | `run_vina.sh` | 30 条 `vina` 调用；盒子中心取自共晶配体质心（无则蛋白质几何中心）；`--exhaustiveness 32 --num_modes 9` | 可执行 |
| 调度脚本 | `run_vina.slurm` | SLURM 头（`--nodes=1 --ntasks-per-node=8 --time=12:00:00 --mem=16G`）；`module load vina` 后调用 `run_vina.sh` | 可执行 |
| 制备脚本 | `prepare_docking.py` | 生成上述全部 PDBQT + manifest + run 脚本；含 SEC/中文路径修复（ASCII 临时目录跑 obabel） | 已运行产生产物 |
| LINCS 重算脚本 | `../02_scripts/step0b_lincs_g1_controls.py` | 用完整 LINCS L1000 重算 G1 阳性对照恢复闸门 | 干跑通过（见 §3） |
| LINCS 对照清单 | `../01_data/lincs/control_sig_ids.txt` | 7 个阳性对照药在 L1000 中的 **112** 个 sig_id | 112 行 |
| LINCS HPC 下载 | `../02_scripts/download_lincs_lvl3_hpc.sh` | HPC 上断点续传完整 Level3 矩阵（~5 GB） | 待 HPC 执行 |

**重新校验命令（在 HPC 或本地均可）：**
```bash
echo "受体 PDBQT 大小（应全部 >0）："; ls -l rec_pdbqt/*.pdbqt | awk '{print $5, $9}'
echo "配体 PDBQT 大小（应全部 >0）："; ls -l lig_pdbqt/*.pdbqt | awk '{print $5, $9}'
echo "PDBQT 扭转树完整性（每个配体应含 ROOT/H PARENT）："; grep -l "ROOT" lig_pdbqt/*.pdbqt
```

---

## 2. 受体 / 配体 / 盒子中心（溯源，供 HPC 复现）

盒子中心取自各结构共晶配体质心（MALT1/IRAK1 无共晶配体或口袋不明确 → `exploratory_box=True`，盒子放大至 30 Å 兼顾探索）：

| 受体 | PDB | 配体（共识 NF-κB 轴） | 盒子中心 (x,y,z) | size (Å) | 探索性? |
|---|---|---|---|---|---|
| JAK3 | 6RU9 | tofa/ruxo/phena/amlex/vorinostat | 124.808, -27.191, 88.508 | 22 | 否 |
| IKBKB | 4KIK | 同上 | 11.267, -6.848, -67.327 | 22 | 否 |
| TBK1 | 4DLE | 同上 | 31.450, -25.030, -12.456 | 22 | 否 |
| MYD88 | 5NFO | 同上 | -24.347, 15.754, -19.622 | 22 | 否 |
| MALT1 | 5CQY | 同上 | 67.689, -63.827, 22.332 | 30 | **是** |
| IRAK1 | 4U2U | 同上 | -29.100, -13.915, 15.876 | 30 | **是** |

> 共晶参照：JAK3_6RU9（tofacitinib 共晶）、IKBKB_4KIK（已报道小分子）、TBK1_4DLE、MYD88_5NFO、IRAK1_4U2U 均取原配体；MALT1_5CQY 与 IRAK1_4U2U 的口袋以结构域几何中心近似，故标注探索性。

---

## 3. Step 0b — LINCS G1 重算（升档路径 a）

**为什么需要：** v3 用的 49,216 签名 LINCS 子集**不含任何 7 个阳性对照药实例** → G1 NOT MET（覆盖度限制，非方法失败）。
**修复：** 在 HPC 用**完整 L1000 矩阵**（GSE92742 Level3，含全部 ~20k 化合物，必然含对照药）重算。

**已验证：** `step0b_lincs_g1_controls.py` 在 49k 子集上干跑，代码路径正确，确认 **0/112 对照命中**（符合预期，证明"子集缺对照"假设成立、重算逻辑无误）。
**对照实例清单：** `01_data/lincs/control_sig_ids.txt` 共 **112** 个 sig_id（lidocaine 33 / memantine 26 / ibudilast 23 / norketamine 19 / gabapentin 4 / paracetamol 4 / minocycline 3）。

**HPC 执行（两步）：**
```bash
# 步骤 1：下载完整 LINCS Level3（HPC 开放网络，~5GB，断点续传）
bash 02_scripts/download_lincs_lvl3_hpc.sh 01_data/lincs
#   -> 产出 01_data/lincs/GSE92742_LVL3.gctx

# 步骤 2：重算 G1
python 02_scripts/step0b_lincs_g1_controls.py \
       --gctx 01_data/lincs/GSE92742_LVL3.gctx \
       --controls 01_data/lincs/control_sig_ids.txt
#   -> 产出 03_results/step0b_summary.json（g1_pass / MWU 单尾 p / 逐药 CS）
```
**判定：** `g1_pass = (MWU 单尾 p<0.05) AND (≥2 对照药 CS<-50)`。若 PASS，则 Step 7/8 重定位由"探索性"升为"闸门验证通过"。

---

## 4. Step 8b — 真实 Vina 对接（升档路径 b）

**HPC 执行：**
```bash
# 方式 A：SLURM 调度（推荐，8 核 / 12h / 16G）
sbatch 06_docking/run_vina.slurm
#   -> 需 HPC 已 `module load vina`（或 conda 环境提供 vina）

# 方式 B：直接 bash（交互/无调度器）
cd 06_docking && bash run_vina.sh
#   -> 产物落 06_docking/results/*.log（每任务最佳打分）+ *.pdbqt（pose）
```

**预期产物（HPC 回传后方可引用）：**
- `06_docking/results/<TARGET>__<LIGAND>.log` — 9 个 pose 的 binding affinity（kcal/mol）。
- `06_docking/results/<TARGET>__<LIGAND>.pdbqt` — 最佳 pose 坐标。
- 汇总脚本（待 HPC 跑后补）：解析 30 个 log 取最佳 pose 打分，生成 `03_results/step8b_vina.csv`。

**真实产物已回传（更新 2026-10-09）：** 30 个 real-box 任务在本地 Vina 1.2.7 重跑完成（与 null 同版本），headline JAK3–amlexanox −8.53 kcal/mol；`decoy_control/analyze_decoy_control.py` 解析 630 个任务并产出 `step8c_decoy_summary.json`；手稿 §3.10/§5/§6 已据真实产物回填，所有 kcal/mol 与 p 值均派生自 `decoy_control/results_*/` 与 `step8c_decoy_summary.json`，无硬编码、无编造。

### 4.1 decoy_control — 特异性对照（已完成）

`06_docking/decoy_control/` 在**本地 Vina 1.2.7**（与 real-box 同版本，避免版本混杂）完成 decoy/random-box 特异性对照：

| 对照类型 | 任务数 | 说明 | 关键结果 |
|---|---|---|---|
| real-box（重跑） | 30 | 6 受体 × 5 配体，40 Å 受体质心盒 | 全部负值；best −8.53 (JAK3–amlexanox)；mean −6.62；10/30 ≤ −7；0 正值 |
| random-box null | 360 | 每受体 12 个 off-target 40 Å 盒（60 盒 × 6 受体） | mean −6.14；min −9.49；6/360 ≤ −8.53 → 精确单尾 p=0.019 |
| decoy-ligand null | 240 | 40 个性质匹配非结合 PDBQT × 6 受体（8/母药 × 5 母药） | amlexanox 组 n=48、mean −6.02、0/48 ≤ −8.53 → p=0.020 |

**解析脚本：** `decoy_control/analyze_decoy_control.py`（纯标准库）→ `decoy_control/step8c_decoy_summary.json`。
**图件：** `04_figures/Fig9_decoy_null.png/.pdf`（双面板 null 分布直方图 + headline 标线 + p 值）。
**诚实结论：** headline 落在两 null 分布尾部（均 p≈0.02，显著但温和）；随机盒最小 −9.49 揭示 40 Å 盒宽容，故对接信号真实但非高置信结合；证据层级仍中档 SCIE（升档硬条件仍非充分满足）。

---

## 5. 诚实边界小结

| 项目 | 状态 | 真实产物 |
|---|---|---|
| Step 0b G1 重算 | ✅ 已完成（HPC 真实执行） | 382/382 对照命中；MWU p=0.0394；0 达 CS<−50 → G1 NOT MET（弱方向信号） |
| Step 8b Vina 对接 | ✅ 已完成（本地 1.2.7 重跑 + decoy 对照） | real-box best −8.53；30/30 负值；randbox null p=0.019；decoy null p=0.020 |
| 手稿 v4 升档 | ✅ 已据真实产物回填 | §3.10/§5/§6 全部更新；无 "decoy pending" 占位 |

**结论：** 两条升档路径已**从"设想"变为"真实执行 + 已回填"**。G1 仍为 NOT MET（弱方向信号）；Step 8b 特异性现已由 decoy 对照支持（p≈0.02，显著但温和），但 40 Å 盒宽容、无实验验证，升档硬条件仍非充分满足。手稿证据层级仍为**稳健中档 SCIE**，不得投发现型 Q1；待实验验证确认结合后方可重评更高档。
