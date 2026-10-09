#!/bin/bash
# ============================================================================
# RUNME.sh — 方案二 升档路径 (a)+(b) HPC 一键执行器
# 依赖：开放网络（下载完整 LINCS L1000 Level3，~5GB）+ vina 二进制
# 用法：在 HPC 节点解压本包后执行  bash RUNME.sh
# 本包已镜像仓库目录结构（00_pipeline/ 01_data/ 02_scripts/ 03_results/ 06_docking/），
# 故 config.py 的相对路径与 run_vina.sh 的相对路径均可正确解析。
# ============================================================================
set -u
BUNDLE=$(cd "$(dirname "$0")" && pwd)
cd "$BUNDLE"
PY=${PYTHON:-python3}
echo "=== 方案二 HPC 升档执行器 (bundle root: $BUNDLE) ==="

# ---------- (a) Step 0b G1 阳性对照恢复闸门重算 ----------
echo
echo "########## [A] Step 0b — LINCS G1 positive-control recovery (full matrix) ##########"
GCTX="$BUNDLE/01_data/lincs/GSE92742_LVL3.gctx"
if [ ! -f "$GCTX" ]; then
  echo "[A.1] 完整 LINCS 矩阵缺失，开始 HPC 下载（断点续传，需开放网络 ~5GB）..."
  bash 02_scripts/download_lincs_lvl3_hpc.sh 01_data/lincs
  GCTX="$BUNDLE/01_data/lincs/GSE92742_LVL3.gctx"
fi
if [ -f "$GCTX" ]; then
  echo "[A.2] 重算 G1（MWU 单尾 + ≥2 药 CS<-50）..."
  $PY 02_scripts/step0b_lincs_g1_controls.py \
      --gctx "$GCTX" \
      --controls 01_data/lincs/control_sig_ids.txt
  echo "[A] 完成 -> 03_results/step0b_summary.json"
  echo "    判定: g1_pass = (MWU 单尾 p<0.05) AND (≥2 对照药 CS<-50)"
else
  echo "[A] 跳过：未提供 LINCS 完整矩阵。可设环境变量 GCTX=/abs/path/to.gctx 后重跑本步。"
fi

# ---------- (b) Step 8b AutoDock Vina 真实对接 ----------
echo
echo "########## [B] Step 8b — AutoDock Vina docking (6 rec x 5 lig = 30 tasks) ##########"
if ! command -v vina >/dev/null 2>&1; then
  echo "[B] 警告: vina 不在 PATH。请先 'module load vina'（或 conda activate 含 vina 的环境），"
  echo "        再单独执行:  bash 06_docking/run_vina.sh   （或 sbatch 06_docking/run_vina.slurm）"
else
  if command -v sbatch >/dev/null 2>&1; then
    echo "[B] 检测到 SLURM -> sbatch 06_docking/run_vina.slurm"
    sbatch 06_docking/run_vina.slurm
  else
    echo "[B] 无调度器 -> bash 06_docking/run_vina.sh"
    bash 06_docking/run_vina.sh
  fi
  echo "[B] 产物 -> 06_docking/results/<TARGET>__<LIGAND>.log （每任务最佳 kcal/mol）"
fi

echo
echo "=== 全部步骤结束。回填手稿 v4 §3.7'/§3.10/闸门表时，仅引用上述真实产物数字，严禁编造。 ==="
