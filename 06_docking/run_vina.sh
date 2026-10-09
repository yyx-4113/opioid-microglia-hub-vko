#!/bin/bash
# AutoDock Vina batch — run on HPC (SCNet/超算). vina must be in PATH.
# Usage: bash run_vina.sh  (after module load / conda activate providing vina)
set -e
DOCK=$(dirname $(readlink -f $0))
OUT="$DOCK/results"
mkdir -p $OUT
EXHAUSTIVE=32   # increase on HPC for thorough sampling

# JAK3_6RU9 x tofacitinib  (JAK3/JAK1 inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/JAK3_6RU9.pdbqt \
      --ligand $DOCK/lig_pdbqt/tofacitinib.pdbqt \
      --center_x 124.808 --center_y -27.191 --center_z 88.508 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/JAK3_6RU9__tofacitinib.pdbqt --log $OUT/JAK3_6RU9__tofacitinib.log

# JAK3_6RU9 x ruxolitinib  (JAK1/2 inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/JAK3_6RU9.pdbqt \
      --ligand $DOCK/lig_pdbqt/ruxolitinib.pdbqt \
      --center_x 124.808 --center_y -27.191 --center_z 88.508 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/JAK3_6RU9__ruxolitinib.pdbqt --log $OUT/JAK3_6RU9__ruxolitinib.log

# JAK3_6RU9 x phenazone  (LINCS reversal hit (CS=-1.45) + CNS-penetrant analgesic)
vina --receptor $DOCK/rec_pdbqt/JAK3_6RU9.pdbqt \
      --ligand $DOCK/lig_pdbqt/phenazone.pdbqt \
      --center_x 124.808 --center_y -27.191 --center_z 88.508 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/JAK3_6RU9__phenazone.pdbqt --log $OUT/JAK3_6RU9__phenazone.log

# JAK3_6RU9 x amlexanox  (TBK1/IKKε inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/JAK3_6RU9.pdbqt \
      --ligand $DOCK/lig_pdbqt/amlexanox.pdbqt \
      --center_x 124.808 --center_y -27.191 --center_z 88.508 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/JAK3_6RU9__amlexanox.pdbqt --log $OUT/JAK3_6RU9__amlexanox.log

# JAK3_6RU9 x vorinostat  (LINCS reversal hit (CS=-1.31); HDAC inhibitor)
vina --receptor $DOCK/rec_pdbqt/JAK3_6RU9.pdbqt \
      --ligand $DOCK/lig_pdbqt/vorinostat.pdbqt \
      --center_x 124.808 --center_y -27.191 --center_z 88.508 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/JAK3_6RU9__vorinostat.pdbqt --log $OUT/JAK3_6RU9__vorinostat.log

# IKBKB_4KIK x tofacitinib  (JAK3/JAK1 inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/IKBKB_4KIK.pdbqt \
      --ligand $DOCK/lig_pdbqt/tofacitinib.pdbqt \
      --center_x 11.267 --center_y -6.848 --center_z -67.327 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IKBKB_4KIK__tofacitinib.pdbqt --log $OUT/IKBKB_4KIK__tofacitinib.log

# IKBKB_4KIK x ruxolitinib  (JAK1/2 inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/IKBKB_4KIK.pdbqt \
      --ligand $DOCK/lig_pdbqt/ruxolitinib.pdbqt \
      --center_x 11.267 --center_y -6.848 --center_z -67.327 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IKBKB_4KIK__ruxolitinib.pdbqt --log $OUT/IKBKB_4KIK__ruxolitinib.log

# IKBKB_4KIK x phenazone  (LINCS reversal hit (CS=-1.45) + CNS-penetrant analgesic)
vina --receptor $DOCK/rec_pdbqt/IKBKB_4KIK.pdbqt \
      --ligand $DOCK/lig_pdbqt/phenazone.pdbqt \
      --center_x 11.267 --center_y -6.848 --center_z -67.327 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IKBKB_4KIK__phenazone.pdbqt --log $OUT/IKBKB_4KIK__phenazone.log

# IKBKB_4KIK x amlexanox  (TBK1/IKKε inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/IKBKB_4KIK.pdbqt \
      --ligand $DOCK/lig_pdbqt/amlexanox.pdbqt \
      --center_x 11.267 --center_y -6.848 --center_z -67.327 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IKBKB_4KIK__amlexanox.pdbqt --log $OUT/IKBKB_4KIK__amlexanox.log

# IKBKB_4KIK x vorinostat  (LINCS reversal hit (CS=-1.31); HDAC inhibitor)
vina --receptor $DOCK/rec_pdbqt/IKBKB_4KIK.pdbqt \
      --ligand $DOCK/lig_pdbqt/vorinostat.pdbqt \
      --center_x 11.267 --center_y -6.848 --center_z -67.327 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IKBKB_4KIK__vorinostat.pdbqt --log $OUT/IKBKB_4KIK__vorinostat.log

# MALT1_5CQY x tofacitinib  (JAK3/JAK1 inhibitor; targets converged NF-kB axis)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/MALT1_5CQY.pdbqt \
      --ligand $DOCK/lig_pdbqt/tofacitinib.pdbqt \
      --center_x 67.689 --center_y -63.827 --center_z 22.332 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MALT1_5CQY__tofacitinib.pdbqt --log $OUT/MALT1_5CQY__tofacitinib.log

# MALT1_5CQY x ruxolitinib  (JAK1/2 inhibitor; targets converged NF-kB axis)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/MALT1_5CQY.pdbqt \
      --ligand $DOCK/lig_pdbqt/ruxolitinib.pdbqt \
      --center_x 67.689 --center_y -63.827 --center_z 22.332 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MALT1_5CQY__ruxolitinib.pdbqt --log $OUT/MALT1_5CQY__ruxolitinib.log

# MALT1_5CQY x phenazone  (LINCS reversal hit (CS=-1.45) + CNS-penetrant analgesic)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/MALT1_5CQY.pdbqt \
      --ligand $DOCK/lig_pdbqt/phenazone.pdbqt \
      --center_x 67.689 --center_y -63.827 --center_z 22.332 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MALT1_5CQY__phenazone.pdbqt --log $OUT/MALT1_5CQY__phenazone.log

# MALT1_5CQY x amlexanox  (TBK1/IKKε inhibitor; targets converged NF-kB axis)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/MALT1_5CQY.pdbqt \
      --ligand $DOCK/lig_pdbqt/amlexanox.pdbqt \
      --center_x 67.689 --center_y -63.827 --center_z 22.332 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MALT1_5CQY__amlexanox.pdbqt --log $OUT/MALT1_5CQY__amlexanox.log

# MALT1_5CQY x vorinostat  (LINCS reversal hit (CS=-1.31); HDAC inhibitor)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/MALT1_5CQY.pdbqt \
      --ligand $DOCK/lig_pdbqt/vorinostat.pdbqt \
      --center_x 67.689 --center_y -63.827 --center_z 22.332 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MALT1_5CQY__vorinostat.pdbqt --log $OUT/MALT1_5CQY__vorinostat.log

# TBK1_4DLE x tofacitinib  (JAK3/JAK1 inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/TBK1_4DLE.pdbqt \
      --ligand $DOCK/lig_pdbqt/tofacitinib.pdbqt \
      --center_x 31.450 --center_y -25.030 --center_z -12.456 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/TBK1_4DLE__tofacitinib.pdbqt --log $OUT/TBK1_4DLE__tofacitinib.log

# TBK1_4DLE x ruxolitinib  (JAK1/2 inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/TBK1_4DLE.pdbqt \
      --ligand $DOCK/lig_pdbqt/ruxolitinib.pdbqt \
      --center_x 31.450 --center_y -25.030 --center_z -12.456 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/TBK1_4DLE__ruxolitinib.pdbqt --log $OUT/TBK1_4DLE__ruxolitinib.log

# TBK1_4DLE x phenazone  (LINCS reversal hit (CS=-1.45) + CNS-penetrant analgesic)
vina --receptor $DOCK/rec_pdbqt/TBK1_4DLE.pdbqt \
      --ligand $DOCK/lig_pdbqt/phenazone.pdbqt \
      --center_x 31.450 --center_y -25.030 --center_z -12.456 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/TBK1_4DLE__phenazone.pdbqt --log $OUT/TBK1_4DLE__phenazone.log

# TBK1_4DLE x amlexanox  (TBK1/IKKε inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/TBK1_4DLE.pdbqt \
      --ligand $DOCK/lig_pdbqt/amlexanox.pdbqt \
      --center_x 31.450 --center_y -25.030 --center_z -12.456 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/TBK1_4DLE__amlexanox.pdbqt --log $OUT/TBK1_4DLE__amlexanox.log

# TBK1_4DLE x vorinostat  (LINCS reversal hit (CS=-1.31); HDAC inhibitor)
vina --receptor $DOCK/rec_pdbqt/TBK1_4DLE.pdbqt \
      --ligand $DOCK/lig_pdbqt/vorinostat.pdbqt \
      --center_x 31.450 --center_y -25.030 --center_z -12.456 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/TBK1_4DLE__vorinostat.pdbqt --log $OUT/TBK1_4DLE__vorinostat.log

# IRAK1_4U2U x tofacitinib  (JAK3/JAK1 inhibitor; targets converged NF-kB axis)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/IRAK1_4U2U.pdbqt \
      --ligand $DOCK/lig_pdbqt/tofacitinib.pdbqt \
      --center_x -29.100 --center_y -13.915 --center_z 15.876 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IRAK1_4U2U__tofacitinib.pdbqt --log $OUT/IRAK1_4U2U__tofacitinib.log

# IRAK1_4U2U x ruxolitinib  (JAK1/2 inhibitor; targets converged NF-kB axis)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/IRAK1_4U2U.pdbqt \
      --ligand $DOCK/lig_pdbqt/ruxolitinib.pdbqt \
      --center_x -29.100 --center_y -13.915 --center_z 15.876 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IRAK1_4U2U__ruxolitinib.pdbqt --log $OUT/IRAK1_4U2U__ruxolitinib.log

# IRAK1_4U2U x phenazone  (LINCS reversal hit (CS=-1.45) + CNS-penetrant analgesic)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/IRAK1_4U2U.pdbqt \
      --ligand $DOCK/lig_pdbqt/phenazone.pdbqt \
      --center_x -29.100 --center_y -13.915 --center_z 15.876 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IRAK1_4U2U__phenazone.pdbqt --log $OUT/IRAK1_4U2U__phenazone.log

# IRAK1_4U2U x amlexanox  (TBK1/IKKε inhibitor; targets converged NF-kB axis)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/IRAK1_4U2U.pdbqt \
      --ligand $DOCK/lig_pdbqt/amlexanox.pdbqt \
      --center_x -29.100 --center_y -13.915 --center_z 15.876 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IRAK1_4U2U__amlexanox.pdbqt --log $OUT/IRAK1_4U2U__amlexanox.log

# IRAK1_4U2U x vorinostat  (LINCS reversal hit (CS=-1.31); HDAC inhibitor)  [EXPLORATORY BOX]
vina --receptor $DOCK/rec_pdbqt/IRAK1_4U2U.pdbqt \
      --ligand $DOCK/lig_pdbqt/vorinostat.pdbqt \
      --center_x -29.100 --center_y -13.915 --center_z 15.876 \
      --size_x 30.0 --size_y 30.0 --size_z 30.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/IRAK1_4U2U__vorinostat.pdbqt --log $OUT/IRAK1_4U2U__vorinostat.log

# MYD88_5NFO x tofacitinib  (JAK3/JAK1 inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/MYD88_5NFO.pdbqt \
      --ligand $DOCK/lig_pdbqt/tofacitinib.pdbqt \
      --center_x -24.347 --center_y 15.754 --center_z -19.622 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MYD88_5NFO__tofacitinib.pdbqt --log $OUT/MYD88_5NFO__tofacitinib.log

# MYD88_5NFO x ruxolitinib  (JAK1/2 inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/MYD88_5NFO.pdbqt \
      --ligand $DOCK/lig_pdbqt/ruxolitinib.pdbqt \
      --center_x -24.347 --center_y 15.754 --center_z -19.622 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MYD88_5NFO__ruxolitinib.pdbqt --log $OUT/MYD88_5NFO__ruxolitinib.log

# MYD88_5NFO x phenazone  (LINCS reversal hit (CS=-1.45) + CNS-penetrant analgesic)
vina --receptor $DOCK/rec_pdbqt/MYD88_5NFO.pdbqt \
      --ligand $DOCK/lig_pdbqt/phenazone.pdbqt \
      --center_x -24.347 --center_y 15.754 --center_z -19.622 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MYD88_5NFO__phenazone.pdbqt --log $OUT/MYD88_5NFO__phenazone.log

# MYD88_5NFO x amlexanox  (TBK1/IKKε inhibitor; targets converged NF-kB axis)
vina --receptor $DOCK/rec_pdbqt/MYD88_5NFO.pdbqt \
      --ligand $DOCK/lig_pdbqt/amlexanox.pdbqt \
      --center_x -24.347 --center_y 15.754 --center_z -19.622 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MYD88_5NFO__amlexanox.pdbqt --log $OUT/MYD88_5NFO__amlexanox.log

# MYD88_5NFO x vorinostat  (LINCS reversal hit (CS=-1.31); HDAC inhibitor)
vina --receptor $DOCK/rec_pdbqt/MYD88_5NFO.pdbqt \
      --ligand $DOCK/lig_pdbqt/vorinostat.pdbqt \
      --center_x -24.347 --center_y 15.754 --center_z -19.622 \
      --size_x 22.0 --size_y 22.0 --size_z 22.0 \
      --exhaustiveness $EXHAUSTIVE --num_modes 9 \
      --out $OUT/MYD88_5NFO__vorinostat.pdbqt --log $OUT/MYD88_5NFO__vorinostat.log
