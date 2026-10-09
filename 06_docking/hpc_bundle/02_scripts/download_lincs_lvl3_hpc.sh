#!/bin/bash
# ============================================================================
# download_lincs_lvl3_hpc.sh
# 在 HPC (SCNet / 超算，开放网络) 上下载完整 LINCS L1000 Level 3 矩阵，
# 供 step0b_lincs_g1_controls.py 抽取 7 个阳性对照药实例、重算 G1。
# 沙箱代理会拦截 GEO 大文件，故须在 HPC 执行本脚本。
# ============================================================================
set -e
OUT_DIR="${1:-01_data/lincs}"
mkdir -p "$OUT_DIR"
cd "$(dirname "$0")/.."   # 回到项目根
OUT_DIR="$(pwd)/$OUT_DIR"

# GSE92742 完整 Level 3 (z-scored landmark, 1.3M 签名 x 978) — 约 5 GB，用断点续传
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x978.gctx.gz"
DEST="$OUT_DIR/GSE92742_LVL3.gctx.gz"

echo "[download] 目标: $DEST"
if [ -f "$DEST" ]; then
  echo "[download] 已存在，尝试断点续传..."
  curl -C - -L --retry 5 --retry-delay 10 -o "$DEST" "$URL"
else
  curl -L --retry 5 --retry-delay 10 -o "$DEST" "$URL"
fi

echo "[download] 解压..."
gunzip -f -k "$DEST"   # 生成 GSE92742_LVL3.gctx

echo "[download] 完成。随后运行："
echo "python 02_scripts/step0b_lincs_g1_controls.py --gctx $OUT_DIR/GSE92742_LVL3.gctx"
