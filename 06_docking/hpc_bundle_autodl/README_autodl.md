# 方案二 升档路径 (a)+(b) — AUTODL 执行包使用说明

本包是**自包含、可移植**的，用于在 **AUTODL**（https://www.autodl.com）云算力实例上
真实执行两条升档路径，把重定位线从"探索性"升为"硬阳性"：

- **(a) Step 0b G1 重算**：下载完整 LINCS L1000 Level3 矩阵（~5GB），抽取 7 个
  已知阳性对照药（minocycline / ibudilast / norketamine / memantine / gabapentin /
  paracetamol / lidocaine）命中的 **112 个 sig_id**，用全矩阵重算 G1 闸门
  （MWU 单尾 p<0.05 且 ≥2 对照药 CS<-50）。
- **(b) Step 8b 真实 Vina 对接**：6 个靶点受体（JAK3_6RU9 / IKBKB_4KIK / MALT1_5CQY /
  TBK1_4DLE / IRAK1_4U2U / MYD88_5NFO） × 5 个配体（tofacitinib / ruxolitinib /
  phenazone / amlexanox / vorinostat）= 30 个对接任务，结果落 `06_docking/results/`。

> 沙箱（WorkBuddy 对话环境）受代理限制**无法**下载 GEO 大文件、且 `pip install vina`
> 因缺 Boost 编译失败，故真实计算必须在本包 + AUTODL 上完成。本包已**不带** ~5GB
> 的 gctx，由 AUTODL 在执行时自下载。

---

## 一、在 AUTODL 租实例

1. 登录 https://www.autodl.com → 控制台 → 租用新实例。
2. **镜像**：选任意**带 conda 的 Linux 镜像**即可（如 `PyTorch` 系列、或 `Miniconda` /
   `Ubuntu + 自带 conda`）。AUTODL 默认镜像多含 conda。
3. **GPU/CPU**：Vina 用 CPU 即可跑完（30 任务，每任务数秒~数分钟）。
   想更快可租 GPU 实例，Vina 1.2+ 会自动利用 GPU；无需改脚本。
4. **磁盘**：系统盘默认约 30GB，下载 5GB LINCS + 产物 <10GB，足够；
   若担心占满，可在租用时挂**数据盘**（如 `/autodl-fs`），并把本包解压到数据盘。
5. 开机后，用实例提供的 **SSH（公网 IP+端口）** 或 **网页终端/Jupyter** 进入。

## 二、把本包传到实例

两种方式任选其一（从你**本机** D 盘那份 `hpc_bundle_autodl.tar.gz` 出发）：

- **方式 A · 网页上传**：在 AUTODL 实例页打开 Jupyter/JupyterLab 或"文件"面板，
  直接把 `hpc_bundle_autodl.tar.gz` 拖进去。
- **方式 B · scp**：AUTODL 给每个实例一个公网地址与端口（形如 `region-1.autodl.com:3xxxx`），
  从本机终端执行（把 `<host>:<port>` 换成实例信息，用户名为 `root`）：
  ```bash
  scp -P <port> hpc_bundle_autodl.tar.gz root@<host>:/root/
  ```

## 三、解包并跑

```bash
# 进入实例工作目录（例如）
cd /root
tar -xzf hpc_bundle_autodl.tar.gz
cd hpc_bundle_autodl

# 1) 一次性装环境（建 conda env + 装 vina + python 依赖，约几分钟）
bash setup_autodl.sh

# 2) 一键跑 (a)+(b)（下载 5GB + step0b + 30 个 Vina 对接，约几十分钟，视实例而定）
bash RUNME_autodl.sh
```

`RUNME_autodl.sh` 会自动：缺 gctx 时先下载（断点续传）→ 跑 step0b 重算 G1 → 跑
`06_docking/run_vina.sh` 完成 30 个对接。中途断网可用 `bash RUNME_autodl.sh` 重跑，
gctx 已存在则跳过下载。

## 四、回传产物（回填手稿 v4）

执行完成后，把以下文件从实例取回（网页下载或 scp），交回 WorkBuddy 会话即可回填：

- `03_results/step0b_summary.json`         —— G1 闸门真实结果（MWU p 值、各对照药 CS、g1_pass）
- `06_docking/results/*.log`               —— 30 个对接任务每任务最佳 kcal/mol（log 末段 `Affinity:`）
- （可选）`06_docking/results/*.pdbqt`     —— 最佳 pose 结构，可用于作图

回填规则（已在手稿 v4 写定，严禁编造）：
- 若 **G1_PASS = True**（MWU 单尾 p<0.05 且 ≥2 对照药 CS<-50）→ G1 闸门由 NOT MET 转 PASS；
- 若 **某真实对接最佳 pose ≤ −7.0 kcal/mol** 且优于 decoy/对照 → 重定位线升为硬阳性；
- 任一成立即触发向**发现型 Q1**（如 JNI / BBI / BJA）的升档重评，
  并改写 v4 的 §3.7' / §3.10 / 闸门表；二者皆未达成则维持稳健中档 SCIE 结论不变。

## 五、目录结构（包内）

```
hpc_bundle_autodl/
├── 00_pipeline/config.py                  # 中央配置（路径相对 ROOT，包内自解析）
├── 01_data/lincs/
│   ├── gene_info.txt.gz                  # landmark 基因映射（step0b 必需）
│   ├── control_sig_ids.txt               # 112 个阳性对照 sig_id
│   └── (运行时生成) GSE92742_LVL3.gctx   # 由 download 脚本在 AUTODL 下载
├── 02_scripts/
│   ├── step0_lincs.py                    # 复用：landmark/签名构建
│   ├── step0b_lincs_g1_controls.py       # (a) G1 重算主程序
│   └── download_lincs_lvl3_hpc.sh        # (a) 完整 LINCS 下载（断点续传）
├── 03_results/
│   ├── step1_activation_up_full.csv      # NP 小胶质激活上调整体签名（step0b 输入）
│   ├── step1_activation_down_full.csv    # NP 小胶质激活下调整体签名
│   └── step0b_summary.json               # (a) 产出（AUTODL 生成）
├── 06_docking/
│   ├── rec_pdbqt/*.pdbqt                 # 6 受体 PDBQT（已生成，全非零）
│   ├── lig_pdbqt/*.pdbqt                 # 5 配体 PDBQT（已生成，含完整扭转树）
│   ├── run_vina.sh                       # 30 任务 Vina 批处理
│   ├── run_vina.slurm                    # （SLURM 参考，AUTODL 不用）
│   ├── manifest.csv                      # 30 任务清单
│   └── results/                          # (b) 产出（AUTODL 生成 *.log）
├── setup_autodl.sh                       # ★ 建环境
├── RUNME_autodl.sh                       # ★ 一键执行 (a)+(b)
└── README_autodl.md
```

> 注：原通用版 `RUNME.sh` / `README_bundle.md` 一并保留作 SLURM/超算参考，
> 在 AUTODL 上请使用 `setup_autodl.sh` + `RUNME_autodl.sh`。
