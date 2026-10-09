#!/usr/bin/env bash
# ============================================================================
# setup_autodl.sh — 在 AUTODL / SeetaCloud 实例上一键建好运行环境
# 实测：Ubuntu 22.04，apt 仓库直接提供 autodock-vina 1.2.3（含二进制，无需 Boost 编译）；
#       conda-forge 国内镜像不可达，故改用 apt + 清华 pip 镜像。
# 流程：
#   1) apt-get 安装 autodock-vina（提供 vina 命令）
#   2) apt-get 安装 python3-pip
#   3) pip 清华镜像安装 h5py/scipy/numpy/pandas（step0b 读 gctx 与 MWU 所需）
# 用法：bash setup_autodl.sh   （只需跑一次）
# ============================================================================
set -e
echo "[setup] 更新 apt 并安装 autodock-vina + python3-pip ..."
apt-get update -y
apt-get install -y autodock-vina python3-pip

echo "[setup] 安装 python 依赖 (h5py/scipy/numpy/pandas, 清华镜像) ..."
pip3 install -i https://pypi.tuna.tsinghua.edu.cn/simple h5py scipy numpy pandas

echo "[setup] 验证安装:"
which vina && (vina --version || true)
python3 -c "import h5py,scipy,numpy,pandas; print('python deps OK')"
echo "[setup] 完成。下一步运行: bash RUNME_autodl.sh"
