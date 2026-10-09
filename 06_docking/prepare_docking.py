# -*- coding: utf-8 -*-
"""
方案二 — Step 8(b): 分子对接输入包制备 (HPC-ready, 沙箱内可跑前处理)
=========================================================================
本脚本在沙箱内完成所有"前处理"（无需 vina 二进制）：
  1. 受体：调用 meeko mk_prepare_receptor.py 将 PDB -> PDBQT，并据共晶配体自动设定对接盒子；
  2. 配体：rdkit 生成 3D 构象 -> meeko 生成 PDBQT；
  3. 生成 HPC 运行脚本 run_vina.sh / run_vina.slurm 与任务清单 manifest.csv。

真实对接 (vina --receptor ... --ligand ... ) 在 HPC (SCNet/超算) 执行，本脚本不伪造打分。

依赖(沙箱内已具备)：rdkit, meeko 0.8.0, openbabel 3.2.1。
"""
import os, sys, csv, subprocess, argparse
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCK = os.path.join(ROOT, "06_docking")
REC_DIR = os.path.join(DOCK, "receptors")
RECPQ = os.path.join(DOCK, "rec_pdbqt")
LIGPQ = os.path.join(DOCK, "lig_pdbqt")
CFG = os.path.join(DOCK, "configs")
RES = os.path.join(DOCK, "results")
for d in (REC_DIR, RECPQ, LIGPQ, CFG, RES):
    os.makedirs(d, exist_ok=True)

PY = sys.executable
OBABEL = os.path.join(os.path.dirname(sys.executable), "obabel.exe")

# 共晶配体 / 目标 (PDB id 已下载到 receptors/)
TARGETS = {
    "JAK3_6RU9":  ("JAK3",  "6RU9"),
    "IKBKB_4KIK": ("IKBKB", "4KIK"),
    "MALT1_5CQY": ("MALT1", "5CQY"),
    "TBK1_4DLE":  ("TBK1",  "4DLE"),
    "IRAK1_4U2U": ("IRAK1", "4U2U"),
    "MYD88_5NFO": ("MYD88", "5NFO"),  # 探索性：TIR 二聚界面，无小配体
}
# 候选重定位化合物（来自 step8 ADMET + LINCS 逆转命中）
LIGANDS = [
    ("tofacitinib",  "C[C@@H]1CCN(C[C@@H]1N(C)C2=NC=NC3=C2C=CN3)C(=O)CC#N", "JAK3/JAK1 inhibitor; targets converged NF-kB axis"),
    ("ruxolitinib",  "C1CCC(C1)[C@@H](CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3", "JAK1/2 inhibitor; targets converged NF-kB axis"),
    ("phenazone",    "CC1=CC(=O)N(N1C)C2=CC=CC=C2",                        "LINCS reversal hit (CS=-1.45) + CNS-penetrant analgesic"),
    ("amlexanox",    "CC(C)C1=CC2=C(C=C1)OC3=NC(=C(C=C3C2=O)C(=O)O)N",     "TBK1/IKKε inhibitor; targets converged NF-kB axis"),
    ("vorinostat",   "C1=CC=C(C=C1)NC(=O)CCCCCCC(=O)NO",                   "LINCS reversal hit (CS=-1.31); HDAC inhibitor"),
]

SKIP_RES = {"HOH", "WAT", "H2O", "NA", "CL", "CA", "MG", "ZN", "MN", "FE", "K", "SO4", "PO4",
            "HPO4", "NO3", "CO3", "DMS", "EDO", "GOL", "PEG", "ACT", "CIT", "FMT", "ACE", "NH4"}

def pdb_ligand_centroid(pdb_path):
    """返回共晶配体重原子质心 (x,y,z, natoms)，若无合适配体返回 (None,...)"""
    from collections import defaultdict as dd
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("HETATM") or line.startswith("ATOM")):
                continue
            res = line[17:20].strip()
            if res in SKIP_RES:
                continue
            if line.startswith("ATOM"):
                continue  # 只看 HETATM 非水/非离子
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            atoms.append((x, y, z))
    if len(atoms) < 4:
        return None, 0
    cx = sum(a[0] for a in atoms) / len(atoms)
    cy = sum(a[1] for a in atoms) / len(atoms)
    cz = sum(a[2] for a in atoms) / len(atoms)
    return (cx, cy, cz), len(atoms)

def protein_centroid(pdb_path):
    xs = ys = zs = n = 0.0
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            xs += x; ys += y; zs += z; n += 1
    if n == 0:
        return (0.0, 0.0, 0.0)
    return (xs / n, ys / n, zs / n)

STANDARD_AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
                "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}

def write_protein_only_pdb(src_pdb, dst_pdb):
    """仅保留标准氨基酸 ATOM（剥离 HETATM 与非标准残基如 SEC/PYL/磷酸化修饰），
    作为刚性受体。openbabel 对非标准残基的键感知会失败（如 JAK3 的 SEC）。"""
    n = 0
    with open(src_pdb) as f, open(dst_pdb, "w") as g:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            res = line[17:20].strip()
            if res not in STANDARD_AA:
                continue
            g.write(line)
            n += 1
    return n

def prep_receptor(name_pdbkey, pdbid):
    pdb = os.path.join(REC_DIR, f"{name_pdbkey}.pdb")
    out_base = os.path.join(RECPQ, name_pdbkey)
    cleaned = os.path.join(RECPQ, name_pdbkey + "_clean.pdb")
    cen, nlig = pdb_ligand_centroid(pdb)
    if cen is None:
        cen = protein_centroid(pdb)
        size = (30.0, 30.0, 30.0)
        exploratory = True
        note = "no co-crystal small ligand -> box at protein centroid (exploratory)"
    else:
        size = (22.0, 22.0, 22.0)
        exploratory = False
        note = f"box centered on co-crystal ligand centroid (n_atoms={nlig})"
    # 仅留标准氨基酸作为刚性受体，剥离水/离子/共晶配体/非标准残基
    natoms = write_protein_only_pdb(pdb, cleaned)
    # openbabel 二进制无法处理中文路径 -> 在 ASCII 临时目录跑，再拷回
    import tempfile, shutil
    tmp = tempfile.mkdtemp(prefix="dock_")
    clean_tmp = os.path.join(tmp, name_pdbkey + "_clean.pdb")
    out_tmp = os.path.join(tmp, name_pdbkey + ".pdbqt")
    shutil.copy(cleaned, clean_tmp)
    cmd = [OBABEL, clean_tmp, "-O", out_tmp, "-h", "-xr", "--partialcharge", "gasteiger"]
    r = subprocess.run(cmd, capture_output=True, encoding="utf-8", errors="ignore")
    ok = os.path.exists(out_tmp) and os.path.getsize(out_tmp) > 0
    if ok:
        shutil.copy(out_tmp, out_base + ".pdbqt")
    print(f"[rec] {name_pdbkey}: {'OK' if ok else 'FAIL'} natoms={natoms} box=({cen[0]:.1f},{cen[1]:.1f},{cen[2]:.1f}) {note}")
    if not ok:
        print("   stderr:", r.stderr[-300:])
    return ok, cen, size, exploratory, note

def prep_ligand(name, smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation
    out = os.path.join(LIGPQ, f"{name}.pdbqt")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[lig] {name}: SMILES parse FAIL"); return False
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    prep.write_pdbqt_file(out)
    ok = os.path.exists(out)
    print(f"[lig] {name}: {'OK' if ok else 'FAIL'} ({len(setups)} setups)")
    return ok

def main():
    rec_status = {}
    for key, (gene, pid) in TARGETS.items():
        ok, cen, size, exp, note = prep_receptor(key, pid)
        rec_status[key] = (ok, cen, size, exp, note, gene)
    lig_ok = {}
    for name, smi, desc in LIGANDS:
        lig_ok[name] = prep_ligand(name, smi)

    # 生成 manifest 与 HPC 运行脚本
    jobs = []
    for key, (gene, pid) in TARGETS.items():
        ok, cen, size, exp, note, g = rec_status[key]
        if not ok:
            continue
        for name, smi, desc in LIGANDS:
            if not lig_ok.get(name):
                continue
            jobs.append((key, g, name, desc, exp))
    # manifest.csv
    with open(os.path.join(DOCK, "manifest.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target_key", "target_gene", "ligand", "ligand_note", "exploratory_box"])
        for j in jobs:
            w.writerow([j[0], j[1], j[2], j[3], j[4]])
    # run_vina.sh
    sh = ["#!/bin/bash", "# AutoDock Vina batch — run on HPC (SCNet/超算). vina must be in PATH.",
          "# Usage: bash run_vina.sh  (after module load / conda activate providing vina)",
          "set -e", "DOCK=$(dirname $(readlink -f $0))", 'OUT="$DOCK/results"', "mkdir -p $OUT",
          'EXHAUSTIVE=32   # increase on HPC for thorough sampling', ""]
    for key, g, name, desc, exp in jobs:
        box = rec_status[key][1]
        sh.append(f"# {key} x {name}  ({desc})" + ("  [EXPLORATORY BOX]" if exp else ""))
        sh.append(f"vina --receptor $DOCK/rec_pdbqt/{key}.pdbqt \\")
        sh.append(f"      --ligand $DOCK/lig_pdbqt/{name}.pdbqt \\")
        sh.append(f"      --center_x {box[0]:.3f} --center_y {box[1]:.3f} --center_z {box[2]:.3f} \\")
        sh.append(f"      --size_x {rec_status[key][2][0]:.1f} --size_y {rec_status[key][2][1]:.1f} --size_z {rec_status[key][2][2]:.1f} \\")
        sh.append(f"      --exhaustiveness $EXHAUSTIVE --num_modes 9 \\")
        sh.append(f"      --out $OUT/{key}__{name}.pdbqt --log $OUT/{key}__{name}.log")
        sh.append("")
    with open(os.path.join(DOCK, "run_vina.sh"), "w") as f:
        f.write("\n".join(sh))
    # slurm
    slurm = ["#!/bin/bash", "#SBATCH --job-name=lincs_vina", "#SBATCH --nodes=1 --ntasks-per-node=8",
             "#SBATCH --time=12:00:00", "#SBATCH --mem=16G", "#SBATCH --output=vina_%A.out",
             "module load vina 2>/dev/null || true", "bash $SLURM_SUBMIT_DIR/run_vina.sh"]
    with open(os.path.join(DOCK, "run_vina.slurm"), "w") as f:
        f.write("\n".join(slurm) + "\n")
    print(f"\n[done] receptors prepared={sum(1 for v in rec_status.values() if v[0])}/{len(rec_status)} "
          f"ligands={sum(lig_ok.values())}/{len(LIGANDS)}  jobs={len(jobs)}")
    print(f"[done] wrote manifest.csv, run_vina.sh, run_vina.slurm to 06_docking/")

if __name__ == "__main__":
    main()
