#!/bin/bash
# run_pipeline.sh — 顺序执行 Step1→2→3→4→5（DESeq2 就绪后重跑，确认 G2/G3 仍 PASS）
# 用法: bash 02_scripts/run_pipeline.sh            # 仅跑 1-5
#       bash 02_scripts/run_pipeline.sh --with-lincs # 额外跑 Step0/7（需 LINCS 矩阵就绪）
set -u
cd "$(dirname "$0")/.."
PY="C:/Users/Administrator/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
export PYTHONPATH="$PWD/00_pipeline:$PYTHONPATH"

run() { echo "========== $1 =========="; "$PY" "02_scripts/$1" || { echo "!! $1 失败，停止"; exit 1; }; }

run step1_deg.py
run step2_enrichment.py
run step3_wgcna.py
run step4_hub.py
run step5_cross_species.py

if [ "${1:-}" = "--with-lincs" ]; then
  if [ -f "01_data/lincs/Level2_GEX_delta.gctx.gz" ]; then
    run step0_lincs.py
  else
    echo "!! LINCS 矩阵缺失，跳过 Step0/7"
  fi
fi
echo "========== PIPELINE DONE =========="
