# -*- coding: utf-8 -*-
"""
方案二 — Step 8(b) decoy / 随机盒 特异性对照 输入制备 (沙箱内可跑前处理)
=========================================================================
产出 (供 HPC 真实 vina 对接, 本脚本不伪造任何打分):
  1. 随机盒清单 random_boxes.csv : 每受体 12 个偏离真实结合位点的 40A 盒
  2. decoy 配体 PDBQT : 5 个活性 -> 各 ~8 个属性匹配 (MW/logP) 的 decoy 分子
  3. vina_cmds.txt + run_decoy_control.sh : 全部 600 条 vina 命令 (360 随机盒 + 240 decoy), 并行跑

真实对照基准 = 03_results/step8b_vina.csv (centroid40A_v2, 30 任务全负值, 最佳 -8.50)
所有对照盒尺寸/对接参数与 v2 真实跑完全一致 (size=40, ex=32, num_modes=9), 仅改变盒位置(decoy 用真实质心盒)或配体(随机盒用真实配体).
"""
import os, sys, csv, random, subprocess
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUNDLE = os.path.join(ROOT, "06_docking", "hpc_bundle_autodl", "06_docking")
REC_PQ = os.path.join(BUNDLE, "rec_pdbqt")
LIG_PQ = os.path.join(BUNDLE, "lig_pdbqt")
OUTDIR = os.path.join(ROOT, "06_docking", "decoy_control")
os.makedirs(OUTDIR, exist_ok=True)
DEC_PQ = os.path.join(OUTDIR, "lig_decoy_pdbqt")
os.makedirs(DEC_PQ, exist_ok=True)

PY = sys.executable
random.seed(20261008)

# 受体 (与 prepare_docking.py 一致) —— 真实 40A 盒 = 受体质心
TARGETS = {
    "JAK3_6RU9":  "JAK3",
    "IKBKB_4KIK": "IKBKB",
    "MALT1_5CQY": "MALT1",
    "TBK1_4DLE":  "TBK1",
    "IRAK1_4U2U": "IRAK1",
    "MYD88_5NFO": "MYD88",
}
# 5 个真实活性配体 (SMILES 取自 prepare_docking.py)
LIGANDS = {
    "tofacitinib": "C[C@@H]1CCN(C[C@@H]1N(C)C2=NC=NC3=C2C=CN3)C(=O)CC#N",
    "ruxolitinib": "C1CCC(C1)[C@@H](CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3",
    "phenazone":   "CC1=CC(=O)N(N1C)C2=CC=CC=C2",
    "amlexanox":   "CC(C)C1=CC2=C(C=C1)OC3=NC(=C(C=C3C2=O)C(=O)O)N",
    "vorinostat":  "C1=CC=C(C=C1)NC(=O)CCCCCCC(=O)NO",
}

# ---------------------------------------------------------------------------
# 1) 受体几何：用 v2 真实 40A 盒心 (产生 step8b_vina.csv 的同一套盒),
#    并解析 rec_pdbqt ATOM 坐标仅用于统计 natoms
# ---------------------------------------------------------------------------
def parse_coords(pdbqt):
    n=0
    with open(pdbqt) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            n+=1
    return n

# v2 真实盒心 (来自 step0b/8b 真实跑的 40A 受体质心盒; 与 step8b_vina.csv 同一套)
V2_CENTERS = {
    "JAK3_6RU9":  (116.5, -21.3,  82.6),
    "IKBKB_4KIK": ( 25.4,  -2.5, -52.4),
    "MALT1_5CQY": ( 67.7, -63.9,  22.3),
    "TBK1_4DLE":  ( 27.6, -23.8,  -9.6),
    "IRAK1_4U2U": (-29.1, -13.9,  15.9),
    "MYD88_5NFO": (-25.5,  14.6, -18.6),
}

REC_INFO={}
for key in TARGETS:
    pq=os.path.join(REC_PQ, key+".pdbqt")
    natoms=parse_coords(pq)
    cen=V2_CENTERS[key]
    REC_INFO[key]=(cen,None)
    print(f"[rec] {key}: v2_box_center=({cen[0]:.1f},{cen[1]:.1f},{cen[2]:.1f}) natoms={natoms}")

# ---------------------------------------------------------------------------
# 2) 随机盒：每受体 12 个, 中心在真实盒心之外 (距真实盒心 24~40A,
#    方向随机). 40A 盒半边=20A, 故 d>20A 保证随机盒不含真实结合位点.
#    不依赖本地包围盒, 避免本地/HPC rec_pdbqt 制备差异导致的坐标错位.
# ---------------------------------------------------------------------------
NB=12
random_boxes=[]
for key in TARGETS:
    cen=REC_INFO[key][0]
    centers=[]
    tries=0
    while len(centers)<NB and tries<20000:
        tries+=1
        # 随机单位向量
        v=(random.gauss(0,1),random.gauss(0,1),random.gauss(0,1))
        norm=(v[0]**2+v[1]**2+v[2]**2)**0.5
        if norm<1e-6: continue
        u=(v[0]/norm,v[1]/norm,v[2]/norm)
        d=random.uniform(24.0,40.0)
        x=cen[0]+u[0]*d; y=cen[1]+u[1]*d; z=cen[2]+u[2]*d
        centers.append((x,y,z))
    if len(centers)<NB:
        print(f"  [warn] {key}: only {len(centers)} off-site boxes found")
    for i,(x,y,z) in enumerate(centers):
        random_boxes.append((key, f"{key}_rand{i:02d}", x,y,z))
print(f"[randbox] total={len(random_boxes)} boxes ({NB}/receptor)")

# ---------------------------------------------------------------------------
# 3) decoy 配体：属性匹配 (MW/logP) 的多样药物样分子, 经 meeko 制备 PDBQT
#    注意：剔除已知 JAK/NF-kB 靶向药与本项目活性/对照药, 仅留预期非结合剂
# ---------------------------------------------------------------------------
DECOY_LIB = [
    ("aspirin","CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen","CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("naproxen","CC(C)Cc1ccc(cc1)C(=O)O"),
    ("diclofenac","O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
    ("warfarin","CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
    ("diazepam","CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
    ("phenytoin","O=C1NC(=O)NC1(c1ccccc1)c1ccccc1"),
    ("baclofen","NCCC(c1ccc(Cl)cc1)C(=O)O"),
    ("caffeine","CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("theophylline","Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
    ("metformin","CN(C)C(=N)NC(=N)N"),
    ("fluoxetine","CNCCCOc1ccc(cc1)C(F)(F)F"),
    ("propranolol","CC(C)NCC(O)COc1cccc2ccccc12"),
    ("atenolol","CC(C)NCC(O)COc1ccc(cc1)NCC(N)=O"),
    ("diphenhydramine","CN(C)CCOC(c1ccccc1)c1ccccc1"),
    ("zidovudine","Cc1c(C)n(c(=O)[nH]1)COC(C)C"),
    ("venlafaxine","CC(NCC(O)CC)C1=CC=C(OCC)C=C1"),
    ("amitriptyline","CN(C)CCC=C1C=CC=CC=C1"),
    ("chlorpromazine","CN(C)CCSC1=CC=C(N(CC)CC)C=C1"),
    ("imipramine","CN(C)CCC1=CC=CC=C1"),
    ("morphine_x","COc1ccc2c(c1)CC(O)C1C2CCC1O"),  # 仅作为化学多样 decoy (非项目活性)
    ("codeine_x","COc1ccc2c(c1OC)CC(O)C1C2CCC1O"),
    ("quinine_x","Cc1cc2[nH]ccc(C(C)CO)c2c1OC"),
    ("ketamine_x","CN1C2CCC1C(N)C2(Cl)Cl"),
    ("loratadine","COc1ccc(cc1)OCCN2CCC(CC2)C(O)c2ccc(Cl)c(c2)C#N"),
    ("fexofenadine","CC(C)(C)c1ccc(cc1)CC(C(=O)O)NCCc2ccc(OCC(C)O)cc2"),
    ("montelukast","COc1ccc(cc1)CC(C(=O)O)NCCc2ccc(cc2)C=Cc3ccc(Cl)cc3"),
    ("omeprazole_x","COc1ccc2[nH]c([nH]c2c1)CS(=O)Cc1ccccc1"),
    ("pantoprazole_x","CS(=O)(=O)c1nc2ccc(CC(=O)O)cc2n1c1ccc(OC)cc1"),
    ("cimetidine","Cc1nc(NC)nc(NCCSCc2ccccc2)n1"),
    ("ranitidine_x","Cc1nc(N)nc(NCCSCc2ccc(CN(C)C)cc2)n1"),
    ("acyclovir_x","Nc1nc2c(ncn2C)NCCO1"),
    ("acyclovir_y","CC(C=O)OC(C)C"),
    ("simvastatin_x","CC(C)C(=O)OC1=C(C)C(C)=C(C)C(=C1)C(=O)CC(C)CC(=O)O"),
    ("atorvastatin_x","CCC(CC)(C)c1ccc(C(C)C(=O)O)cc1C(=O)NCC(O)CC(O)CC(O)CC(O)CC(O)"),
    ("amlodipine_x","CCOC(=C1)CNC1C(=O)NC(=Cc2cc(Cl)ccc2Cl)C(=O)OC2CC(OC)CC(C)O2"),
    ("verapamil_x","COc1ccc(cc1)CC(NC(=O)C(C)COC(=O)c1ccc(OC)c(C)c1)CC(C)N(C)C"),
    ("methotrexate_x","CN(Cc1cnc2c(N)ncnc12)C(=O)N[C@@H](CCC(=O)O)C(=O)O"),
    ("prednisolone_x","CC12CC(O)CC1CC1C2C(=O)CC(O)C1C(=O)C(O)CCC1"),
    ("testosterone_x","CC12CCC3C(C1CCC2O)CCC3C1=CC(=O)CC(C1)C"),
    ("estradiol_x","Cc1cc2c(cc1O)CCC1C2CC(O)CC1"),
    ("caffeine_y","CN1C=NC2=C1C(=O)N(C)C(=O)N2C"),
    ("theobromine","Cn1cnc2c1c(=O)[nH]c(=O)n2C"),
    ("nicotine","CN1CCCC1c1cccnc1"),
    ("amphetamine","CC(N)CCc1ccccc1"),
    ("ephedrine","CC(NC)C(O)Cc1ccccc1"),
]
# 剔除可能与本项目活性/对照冲突的命名后缀 _x/_y 仅作去重, 不影响分子本身

def prep_decoy_pdbqt(name, smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation
    out=os.path.join(DEC_PQ, f"{name}.pdbqt")
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol=Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        prep=MoleculePreparation()
        setups=prep.prepare(mol)
        prep.write_pdbqt_file(out)
    except Exception as e:
        print(f"  [decoy] {name}: prep FAIL ({e})")
        return None
    if os.path.exists(out) and os.path.getsize(out)>0:
        return out
    return None

from rdkit import Chem
from rdkit.Chem import Descriptors
active_props={}
for name,smi in LIGANDS.items():
    m=Chem.MolFromSmiles(smi); m=Chem.AddHs(m)
    active_props[name]=(Descriptors.MolWt(m), Descriptors.MolLogP(m))

# 计算每个 decoy 库分子的 MW/logP, 为每个活性选 8 个最接近 (MW+logP 归一距离)
lib_props={}
for dname,dsmi in DECOY_LIB:
    m=Chem.MolFromSmiles(dsmi)
    if m is None: continue
    m=Chem.AddHs(m)
    lib_props[dname]=(dsmi, Descriptors.MolWt(m), Descriptors.MolLogP(m))

selected_decoys=[]  # (decoy_name, parent_active, smiles)
for aname,(amw,alp) in active_props.items():
    scored=[]
    for dname,(dsmi,dmw,dlp) in lib_props.items():
        if dname==aname: continue
        dmw_n=(dmw-amw)/100.0
        dlp_n=(dlp-alp)/2.0
        dist=(dmw_n**2+dlp_n**2)**0.5
        scored.append((dist,dname,dsmi))
    scored.sort()
    chosen=scored[:8]
    for dist,dname,dsmi in chosen:
        decoy_name=f"decoy_{aname}_{dname}"
        selected_decoys.append((decoy_name,aname,dsmi))
        print(f"[decoy] {aname} <- {dname} (MW={lib_props[dname][1]:.0f},logP={lib_props[dname][2]:.2f})")

# 制备 decoy PDBQT
prep_count=0
for decoy_name,aname,dsmi in selected_decoys:
    if prep_decoy_pdbqt(decoy_name, dsmi):
        prep_count+=1
print(f"[decoy] prepared {prep_count}/{len(selected_decoys)} PDBQTs -> {DEC_PQ}")

# ---------------------------------------------------------------------------
# 4) 写对照输入文件
# ---------------------------------------------------------------------------
with open(os.path.join(OUTDIR,"random_boxes.csv"),"w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["receptor_key","box_id","center_x","center_y","center_z","size"])
    for key,bid,x,y,z in random_boxes:
        w.writerow([key,bid,f"{x:.3f}",f"{y:.3f}",f"{z:.3f}",40.0])

with open(os.path.join(OUTDIR,"decoy_manifest.csv"),"w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["decoy_name","parent_active","smiles"])
    for decoy_name,aname,dsmi in selected_decoys:
        w.writerow([decoy_name,aname,dsmi])

# ---------------------------------------------------------------------------
# 5) 生成 HPC vina 命令 (全部单行, 供 xargs -P 并行)
# ---------------------------------------------------------------------------
SIZE=40.0
EX=32
# 注意: 路径用环境变量 ($REC_PQ/$LIG_PQ/$DEC_PQ/$OUTDIR), 由 HPC 端 run_decoy_control.sh 定义
cmds=[]
# (a) 随机盒对照: 5 真实配体 x 6 受体 x 12 盒
for key,bid,x,y,z in random_boxes:
    for lig in LIGANDS:
        out=f"$OUTDIR/results_randbox/{bid}__{lig}"
        cmds.append(
            f"vina --receptor $REC_PQ/{key}.pdbqt --ligand $LIG_PQ/{lig}.pdbqt "
            f"--center_x {x:.3f} --center_y {y:.3f} --center_z {z:.3f} "
            f"--size_x {SIZE} --size_y {SIZE} --size_z {SIZE} "
            f"--exhaustiveness {EX} --num_modes 9 --out {out}.pdbqt"
        )
# (b) decoy 对照: 每 decoy x 6 受体 (真实 40A 质心盒)
for decoy_name,aname,dsmi in selected_decoys:
    for key,(cen,bbox) in REC_INFO.items():
        out=f"$OUTDIR/results_decoy/{decoy_name}__{key}"
        cmds.append(
            f"vina --receptor $REC_PQ/{key}.pdbqt --ligand $DEC_PQ/{decoy_name}.pdbqt "
            f"--center_x {cen[0]:.3f} --center_y {cen[1]:.3f} --center_z {cen[2]:.3f} "
            f"--size_x {SIZE} --size_y {SIZE} --size_z {SIZE} "
            f"--exhaustiveness {EX} --num_modes 9 --out {out}.pdbqt"
        )
# (c) 真实盒基线 (同版本 Vina 1.2.7 对照真实跑, 与 decoy/随机盒同版本一致): 5 真实配体 x 6 受体 (v2 真实 40A 盒)
for lig in LIGANDS:
    for key,(cen,bbox) in REC_INFO.items():
        out=f"$OUTDIR/results_realbox/{key}__{lig}"
        cmds.append(
            f"vina --receptor $REC_PQ/{key}.pdbqt --ligand $LIG_PQ/{lig}.pdbqt "
            f"--center_x {cen[0]:.3f} --center_y {cen[1]:.3f} --center_z {cen[2]:.3f} "
            f"--size_x {SIZE} --size_y {SIZE} --size_z {SIZE} "
            f"--exhaustiveness {EX} --num_modes 9 --out {out}.pdbqt"
        )

with open(os.path.join(OUTDIR,"vina_cmds.txt"),"w") as f:
    f.write("\n".join(cmds)+"\n")
with open(os.path.join(OUTDIR,"run_decoy_control.sh"),"w") as f:
    f.write("#!/bin/bash\n")
    f.write("# 在 HPC 上定义真实路径后运行; 上传 decoy_control/ 到 $DOCK/decoy_control/\n")
    f.write("set -e\n")
    f.write("DOCK_ROOT=${1:-/root/opiate_dock/hpc_bundle_autodl/06_docking}\n")
    f.write("export REC_PQ=\"$DOCK_ROOT/rec_pdbqt\"\n")
    f.write("export LIG_PQ=\"$DOCK_ROOT/lig_pdbqt\"\n")
    f.write("export DEC_PQ=\"$DOCK_ROOT/decoy_control/lig_decoy_pdbqt\"\n")
    f.write("export OUTDIR=\"$DOCK_ROOT/decoy_control\"\n")
    f.write("mkdir -p \"$OUTDIR/results_randbox\" \"$OUTDIR/results_decoy\"\n")
    f.write("cat \"$OUTDIR/vina_cmds.txt\" | xargs -P 24 -I{} bash -c '{}'\n")
print(f"[cmds] wrote {len(cmds)} vina commands (randbox={len(random_boxes)*5}, decoy={len(selected_decoys)*len(TARGETS)})")
print(f"[done] outputs in {OUTDIR}")
