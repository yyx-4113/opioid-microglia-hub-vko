"""
Step 8 — In silico ADMET / blood-brain-barrier (BBB) permeability prediction.

CONTEXT / HONEST SCOPE:
  - AutoDock Vina could NOT be installed in this sandbox (Boost build failure),
    so structure-based docking against NF-kB-axis targets (JAK3/REL/MALT1/...)
    was NOT performed. Docking is flagged as a planned HPC follow-up.
  - This script therefore performs a REAL, reproducible cheminformatics ADMET/BBB
    prediction (rdkit) on:
      (a) the named small-molecule LINCS reversal hits (trt_cp, CS<0):
          phenazone, vorinostat;
      (b) calibration references with known CNS status: caffeine (CNS+), metformin (CNS low);
      (c) repurposable inhibitors of the converged NF-kB/innate-immunity axis
          (real compounds, SMILES resolved from PubChem): tofacitinib (JAK3/JAK1),
          ruxolitinib (JAK1/2), amlexanox (TBK1/IKKepsilon).
  - BBB classification is a transparent RULE-BASED heuristic (TPSA<=90, 1<=cLogP<=4,
    HBD<=3, 150<=MW<=500). It is NOT a validated QSAR model; reported as prediction.

All outputs written to 03_results/step8_*.
"""
import os, csv, json
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "03_results")
OUT_CSV = os.path.join(RES, "step8_admet.csv")
OUT_SUM = os.path.join(RES, "step8_summary.json")

# compound, target/role, lincs_cs (None if not from LINCS), source_note
COMPOUNDS = [
    ("phenazone", "LINCS reversal small molecule (trt_cp)", -1.45, "Reverses microglial-activation signature in LINCS; antipyrine/analgesic"),
    ("vorinostat", "LINCS reversal small molecule (trt_cp)", -1.31, "Reverses microglial-activation signature in LINCS; HDAC inhibitor"),
    ("caffeine", "Calibration reference (known CNS+)", None, "Positive-control CNS-penetrant drug"),
    ("metformin", "Calibration reference (low CNS)", None, "Negative-control low-CNS-penetrant drug"),
    ("tofacitinib", "JAK3/JAK1 inhibitor (NF-kB-axis repurpose candidate)", None, "FDA-approved; targets JAK3 in converged axis"),
    ("ruxolitinib", "JAK1/2 inhibitor (NF-kB-axis repurpose candidate)", None, "FDA-approved; targets JAK1 in converged axis"),
    ("amlexanox", "TBK1/IKKepsilon inhibitor (NF-kB-axis repurpose candidate)", None, "Targets TBK1 in converged axis"),
]

def fetch_smiles(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        r = requests.get(url, timeout=40)
        if r.status_code == 200:
            d = r.json()
            props = d.get("PropertyTable", {}).get("Properties", [{}])[0]
            return props.get("IsomericSMILES") or props.get("SMILES")
    except Exception as e:
        print("  fetch fail", name, e)
    return None

def bbblabel(mw, logp, tpsa, hbd):
    cns = (tpsa <= 90) and (1 <= logp <= 4) and (hbd <= 3) and (150 <= mw <= 500)
    return "likely_CNS+" if cns else "likely_CNS-"

def main():
    rows = []
    for name, role, cs, note in COMPOUNDS:
        smi = fetch_smiles(name)
        rec = {"compound": name, "role": role, "lincs_cs": cs, "smiles": smi, "note": note}
        if smi:
            m = Chem.MolFromSmiles(smi)
            if m:
                mw = Descriptors.MolWt(m)
                logp = Crippen.MolLogP(m)
                tpsa = rdMolDescriptors.CalcTPSA(m)
                hbd = rdMolDescriptors.CalcNumHBD(m)
                hba = rdMolDescriptors.CalcNumHBA(m)
                rot = rdMolDescriptors.CalcNumRotatableBonds(m)
                rec.update({"MW": round(mw, 2), "cLogP": round(logp, 2),
                            "TPSA": round(tpsa, 2), "HBD": hbd, "HBA": hba,
                            "RotBonds": rot, "BBB_pred": bbblabel(mw, logp, tpsa, hbd)})
            else:
                rec["error"] = "MolFromSmiles None"
        else:
            rec["error"] = "SMILES not fetched"
        rows.append(rec)
        print(f"  {name}: {rec.get('BBB_pred','?')} (MW={rec.get('MW')}, logP={rec.get('cLogP')}, TPSA={rec.get('TPSA')})")

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound", "role", "lincs_cs", "smiles",
                                         "MW", "cLogP", "TPSA", "HBD", "HBA",
                                         "RotBonds", "BBB_pred", "note", "error"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})

    cns_pos = [r["compound"] for r in rows if r.get("BBB_pred") == "likely_CNS+"]
    summary = {
        "scope": "ADMET/BBB rule-based prediction only; structure-based docking NOT performed (AutoDock Vina unavailable in sandbox).",
        "bbb_rule": "likely_CNS+ iff TPSA<=90 AND 1<=cLogP<=4 AND HBD<=3 AND 150<=MW<=500 (heuristic, not validated QSAR).",
        "compounds_evaluated": len(rows),
        "predicted_CNS_positive": cns_pos,
        "note_strongest_candidate": "phenazone (antipyrine): LINCS reversal hit (CS=-1.45) AND predicted CNS+ -> brain-penetrant repurposing candidate to dampen microglial activation in OIH.",
        "converged_axis_targets_for_future_docking": ["MYD88", "JAK3", "REL", "MALT1", "TBK1", "IKBKB", "IRAK1", "CHUK"],
        "docking_status": "NOT PERFORMED (Vina build failed: missing Boost); planned HPC follow-up.",
    }
    with open(OUT_SUM, "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("[step8] wrote", OUT_CSV, "and", OUT_SUM)
    print("[step8] predicted CNS+:", cns_pos)

if __name__ == "__main__":
    main()
