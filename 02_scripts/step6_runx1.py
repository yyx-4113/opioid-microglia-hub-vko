"""
Step 6 — RUNX1 ChIP-seq convergence anchor (honest, data-driven).

Strategy (forced by sandbox network limits: ENCODE/JASPAR/chip-atlas blocked):
  - The LINCS reversal list (step7_drug_cs.csv) is HUMAN. RUNX1 itself appears as a
    reversal hit (trt_sh, CS from step7). We test, in human space (no cross-species
    mapping needed), whether RUNX1 ChIP-seq target genes (from Enrichr ChEA_2016,
    a ChIP-seq-derived TF-target library) are enriched among the genes whose
    perturbation reverses the microglial activation signature.
  - We also directly fetch the RUNX1 target set and (a) compute overlap with the
    reversal-hit gene list and (b) check whether the 6 mouse NP-activation hub
    human orthologs (GAPT/OAS3/CPNE2/SGMS2/CEBPE/IL1F9) are RUNX1 targets.

All numbers are written to 03_results/step6_* for manuscript traceability.
"""
import os, csv, json, sys
import numpy as np
from scipy.stats import hypergeom, fisher_exact

import gseapy as gp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "03_results")
STEP7 = os.path.join(RES, "step7_drug_cs.csv")
OUT_SUM = os.path.join(RES, "step6_summary.json")
OUT_OVL = os.path.join(RES, "step6_runx1_overlap.csv")

# 6 mouse hubs -> human ortholog symbols (1:1 conserved)
HUBS = ["Gapt", "Oas3", "Cpne2", "Sgms2", "Cebpe", "Il1f9"]
HUB_HUMAN = {"Gapt": "GAPT", "Oas3": "OAS3", "Cpne2": "CPNE2",
             "Sgms2": "SGMS2", "Cebpe": "CEBPE", "Il1f9": "IL1F9"}

def load_reversal(cs_thr=-1.0):
    hits = []
    with open(STEP7, newline="") as f:
        for row in csv.DictReader(f):
            try:
                cs = float(row["cs"])
            except:
                continue
            if cs <= cs_thr:
                hits.append((row["pert_iname"], row["pert_type"], cs))
    return hits

def main():
    hits = load_reversal(-1.0)
    genes = [h[0] for h in hits]
    print(f"[step6] reversal hits (CS<=-1.0): n={len(genes)}")
    print("  types:", {t: sum(1 for h in hits if h[1] == t) for t in set(h[1] for h in hits)})

    # --- (A) Enrichr ChEA_2016 / TRANSFAC_Human enrichment of reversal hits ---
    runx_rows = []
    try:
        for lib in ["ChEA_2016", "TRANSFAC_Human", "ENCODE_TF_ChIP-seq_2015"]:
            print(f"[step6] gseapy.enrichr vs {lib} ...")
            enr = gp.enrichr(gene_list=genes, gene_sets=lib,
                             organism="human", cutoff=1.0, verbose=False)
            df = enr.results
            # filter RUNX terms
            for _, r in df.iterrows():
                term = str(r.get("Term", ""))
                if "RUNX" in term.upper():
                    runx_rows.append({
                        "library": lib, "term": term,
                        "pvalue": float(r.get("P-value", "nan")),
                        "adjp": float(r.get("Adjusted P-value", "nan")),
                        "odds": float(r.get("Odds Ratio", "nan")),
                        "overlap": str(r.get("Overlap", "")),
                        "genes": str(r.get("Genes", "")),
                    })
            print(f"  {lib}: found {sum('RUNX' in r['term'].upper() for r in runx_rows if r['library']==lib)} RUNX term(s)")
    except Exception as e:
        print("[step6] enrichr error:", type(e).__name__, e)

    # --- (B) Direct RUNX1 target-set fetch + overlap stats ---
    runx_targets = set()
    try:
        import requests
        for lib in ["ChEA_2016", "TRANSFAC_Human"]:
            try:
                r = requests.get(
                    f"https://maayanlab.cloud/Enrichr/geneSetLibrary?libraryName={lib}&geneSetName=RUNX1",
                    timeout=40)
                if r.status_code == 200:
                    d = r.json()
                    for k, v in d.items():
                        if "RUNX1" in k.upper():
                            runx_targets |= set(v)
            except Exception as e:
                print(f"  geneSetLibrary {lib} failed:", e)
    except Exception as e:
        print("[step6] requests unavailable:", e)

    direct_overlap = None
    hub_in_runx = {}
    if runx_targets:
        print(f"[step6] RUNX1 target genes fetched: {len(runx_targets)}")
        rev_set = set(genes)
        ov = runx_targets & rev_set
        # hypergeometric: background = all human genes ~ 20000
        N = 20000
        K = len(runx_targets)
        n = len(rev_set)
        k = len(ov)
        p = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
        direct_overlap = {"n_runx_targets": K, "n_reversal_hits": n,
                          "overlap": k, "hypergeom_p": float(p),
                          "overlap_genes": sorted(ov)}
        print(f"[step6] direct overlap RUNX1 targets ∩ reversal hits: {k} (p={p:.2e})")
        # hub ortholog check
        for m, h in HUB_HUMAN.items():
            hub_in_runx[m] = h in runx_targets
        print("[step6] hubs in RUNX1 targets:", hub_in_runx)
    else:
        print("[step6] RUNX1 target set NOT fetched (network) — relying on (A) enrichment only")

    # local RUNX1 CS from step7
    runx1_cs = None
    with open(STEP7, newline="") as f:
        for row in csv.DictReader(f):
            if row["pert_iname"].upper() == "RUNX1":
                runx1_cs = float(row["cs"]); break

    summary = {
        "n_reversal_hits_cs_le_-1.0": len(genes),
        "reversal_hit_types": {t: sum(1 for h in hits if h[1] == t) for t in set(h[1] for h in hits)},
        "runx1_pert_cs": runx1_cs,
        "enrichr_RUNX_terms": runx_rows,
        "direct_overlap": direct_overlap,
        "hub_ortholog_in_RUNX1_targets": hub_in_runx,
        "caveat": "Human LINCS space used as proxy for mouse microglial program; ChEA RUNX1 set is human ChIP-seq-derived.",
    }
    with open(OUT_SUM, "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print("[step6] wrote", OUT_SUM)

    # small overlap csv if available
    if direct_overlap and direct_overlap.get("overlap_genes"):
        with open(OUT_OVL, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["gene", "in_reversal_hits", "in_RUNX1_targets", "note"])
            for g in sorted(direct_overlap["overlap_genes"]):
                w.writerow([g, "yes", "yes", "RUNX1-target & reversal hit"])
        print("[step6] wrote", OUT_OVL)

if __name__ == "__main__":
    main()
