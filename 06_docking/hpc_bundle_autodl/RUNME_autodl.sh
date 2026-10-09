#!/usr/bin/env bash
# ============================================================================
# RUNME_autodl.sh — 方案二 升档路径 (a)+(b) 在 AUTODL 单机一键执行
# 前置：已运行 bash setup_autodl.sh 建好 opiate-dock 环境
# 流程：
#   (a) Step 0b G1 重算：下载完整 LINCS Level3 (~5GB, 断点续传) + 抽 112 对照 sig 重算 G1
#   (b) Step 8b Vina 对接：6 受体 x 5 配体 = 30 任务，结果落 06_docking/results/
# 产物回填：03_results/step0b_summary.json + 06_docking/results/*.log
# 严禁编造：仅引用此处真实产出数字回填手稿 v4。
# ============================================================================
set -u
BUNDLE=$(cd "$(dirname "$0")" && pwd)
cd "$BUNDLE"
ENV=opiate-dock
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda env list | grep -q "^$ENV "; then
  echo "[RUNME] 环境 $ENV 不存在，请先运行: bash setup_autodl.sh"
  exit 1
fi
conda activate $ENV
echo "=== 方案二 AUTODL 升档执行器 (bundle: $BUNDLE) ==="
echo "python: $(which python)   vina: $(which vina)"

# ---------- (a) Step 0b G1 阳性对照恢复闸门重算 ----------
echo
echo "########## [A] Step 0b — LINCS G1 positive-control recovery (full matrix) ##########"
GCTX="$BUNDLE/01_data/lincs/GSE92742_LVL3.gctx"
if [ ! -f "$GCTX" ]; then
  echo "[A.1] 完整 LINCS 矩阵缺失，开始下载（断点续传，需开放网络 ~5GB）..."
  bash 02_scripts/download_lincs_lvl3_hpc.sh 01_data/lincs
  GCTX="$BUNDLE/01_data/lincs/GSE92742_LVL3.gctx"
fi
if [ -f "$GCTX" ]; then
  echo "[A.2] 重算 G1（MWU 单尾 p<0.05 且 ≥2 对照药 CS<-50）..."
  python 02_scripts/step0b_lincs_g1_controls.py \
      --gctx "$GCTX" \
      --controls 01_data/lincs/control_sig_ids.txt
  echo "[A] 完成 -> 03_results/step0b_summary.json"
  echo "    判定: g1_pass = (MWU 单尾 p<0.05) AND (≥2 对照药 CS<-50)"
else
  echo "[A] 跳过：LINCS 完整矩阵仍缺失。可设 GCTX=/abs/path/to.gctx 后重跑本步。"
fi

# ---------- (b) Step 8b AutoDock Vina 真实对接 ----------
echo
echo "########## [B] Step 8b — AutoDock Vina docking (6 rec x 5 lig = 30 tasks) ##########"
echo "[B] 运行 run_vina.sh (EXHAUSTIVE=32, CPU)。如租了 GPU 实例，可改 run_vina.sh 中 EXHAUSTIVE 或启用 vina GPU。"
bash 06_docking/run_vina.sh
echo "[B] 产物 -> 06_docking/results/<TARGET>__<LIGAND>.log （每任务最佳 kcal/mol 在 log 末段）"

echo
echo "=== 全部步骤结束。把 03_results/step0b_summary.json 与 06_docking/results/ 回传即可回填 v4 §3.7'/§3.10/闸门表。 ==="
