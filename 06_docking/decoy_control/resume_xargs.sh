#!/bin/bash
# RESUME runner (proven xargs path, grep bug FIXED).
# - Skip logic now uses awk to extract --out (no grep option-parsing pitfall).
# - xargs -P 8 bash -c '{}' expands $REC_PQ/$LIG_PQ/$DEC_PQ/$OUTDIR naturally.
# - Any task whose --out already exists and is >= 5000 bytes is skipped, so
#   the sweep is resumable and never re-runs completed work.
set -u
ROOT="D:/2026.9/极速交付9月会员日优惠套路/05_多组学+虚拟敲除药物发现/方案二_阿片耐受与痛觉过敏小胶质枢纽基因虚拟敲除"
DC="$ROOT/06_docking/decoy_control"
RESUME="$DC/vina_cmds_resume.txt"
export REC_PQ="$ROOT/06_docking/hpc_bundle_autodl/06_docking/rec_pdbqt"
export LIG_PQ="$ROOT/06_docking/hpc_bundle_autodl/06_docking/lig_pdbqt"
export DEC_PQ="$DC/lig_decoy_pdbqt"
export OUTDIR="$DC"

# --- build resume list: skip completed outputs -------------------------
> "$RESUME"
while IFS= read -r cmd; do
  [ -z "$cmd" ] && continue
  out=$(printf '%s\n' "$cmd" | awk '{for(i=1;i<=NF;i++) if($i=="--out") print $(i+1)}')
  out_exp="${out//\$OUTDIR/$OUTDIR}"   # expand $OUTDIR for the existence check
  if [ -f "$out_exp" ] && [ "$(stat -c%s "$out_exp" 2>/dev/null || echo 0)" -ge 5000 ]; then
    : # already done -> skip
  else
    printf '%s\n' "$cmd" >> "$RESUME"
  fi
done < "$DC/vina_cmds_abs.txt"

N=$(wc -l < "$RESUME")
echo "RESUME_BUILD done: $(date)  remaining_tasks=$N  (complete skipped = $((630 - N)))"
if [ "$N" -eq 0 ]; then
  echo "NOTHING TO RUN — sweep complete."
  exit 0
fi

cat "$RESUME" | xargs -P 8 -I{} bash -c '{}' && true
echo "DONE_RESUME_XARGS $(date)"
