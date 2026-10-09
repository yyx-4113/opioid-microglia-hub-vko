#!/bin/bash
# 在 HPC 上定义真实路径后运行; 上传 decoy_control/ 到 $DOCK/decoy_control/
set -e
DOCK_ROOT=${1:-/root/opiate_dock/hpc_bundle_autodl/06_docking}
export REC_PQ="$DOCK_ROOT/rec_pdbqt"
export LIG_PQ="$DOCK_ROOT/lig_pdbqt"
export DEC_PQ="$DOCK_ROOT/decoy_control/lig_decoy_pdbqt"
export OUTDIR="$DOCK_ROOT/decoy_control"
mkdir -p "$OUTDIR/results_randbox" "$OUTDIR/results_decoy"
cat "$OUTDIR/vina_cmds.txt" | xargs -P 24 -I{} bash -c '{}'
