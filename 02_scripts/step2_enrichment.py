# -*- coding: utf-8 -*-
"""
Step 2 — 阿片特异性签名通路富集（ORA, Enrichr）
输入：03_results/step1_opioid_specific_up.csv / _down.csv
输出：03_results/step2_*_enrich.csv + step2_summary.json
注：Enrichr 本体为人的基因集；小鼠 symbol 与人对Ortholog多数同名可直接映射，
    命中项按"符号重叠"解释，作为探索性通路注释（与 v2 §4.4 探索层一致）。
"""
import os, json, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "00_pipeline"))
import config as C

OUT = C.RESULTS
os.makedirs(OUT, exist_ok=True)

LIBS = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
    "WikiPathway_2023_Human",
]

def run_enrich(genes, tag):
    import gseapy
    rows = []
    for lib in LIBS:
        try:
            res = gseapy.enrichr(gene_list=genes, gene_sets=lib,
                                 outdir=None, cutoff=0.05, verbose=False)
            df = res.res2d if hasattr(res, "res2d") else res
            df = df.copy(); df["Library"] = lib
            rows.append(df)
            print(f"[step2:{tag}] {lib}: {len(df)} terms @FDR<0.05")
        except Exception as e:
            print(f"[step2:{tag}] {lib} FAIL: {e}")
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    # 保留关键列
    cols = [c for c in ["Term","Overlap","P-value","Adjusted P-value","Genes","Library","Odds Ratio","Combined Score"] if c in out.columns]
    out = out[cols] if cols else out
    out = out.sort_values("Adjusted P-value" if "Adjusted P-value" in out.columns else "P-value")
    out.to_csv(os.path.join(OUT, f"step2_{tag}_enrich.csv"), index=False)
    return out

def top_terms(df, lib, n=8):
    if df is None or len(df) == 0: return []
    sub = df[df["Library"] == lib] if "Library" in df.columns else df
    if len(sub) == 0: return []
    out = sub.head(n)[["Term", "Adjusted P-value"]].copy()
    out["Adjusted P-value"] = out["Adjusted P-value"].map(lambda v: round(float(v), 10))
    return out.to_dict("records")

def main():
    # 用完整 NP 激活签名（重构后主签名）做富集，才有生物学解释力
    uf = os.path.join(OUT, "step1_activation_up_full.csv")
    df = os.path.join(OUT, "step1_activation_down_full.csv")
    up = pd.read_csv(uf if os.path.exists(uf) else os.path.join(OUT,"step1_activation_up.csv"))["gene"].dropna().astype(str).tolist()
    dn = pd.read_csv(df if os.path.exists(df) else os.path.join(OUT,"step1_activation_down.csv"))["gene"].dropna().astype(str).tolist()
    print(f"[step2] up={len(up)} down={len(dn)}")
    up_df = run_enrich(up, "up")
    dn_df = run_enrich(dn, "down")
    summary = {
        "n_up": len(up), "n_down": len(dn),
        "up_sig_terms": int(len(up_df)), "down_sig_terms": int(len(dn_df)),
        "up_GO_BP_top": top_terms(up_df, "GO_Biological_Process_2023"),
        "up_KEGG_top": top_terms(up_df, "KEGG_2021_Human"),
        "down_GO_BP_top": top_terms(dn_df, "GO_Biological_Process_2023"),
        "down_KEGG_top": top_terms(dn_df, "KEGG_2021_Human"),
        "up_top10": (up_df.head(10)[["Term","Library","Adjusted P-value"]].to_dict("records") if len(up_df) else []),
        "down_top10": (dn_df.head(10)[["Term","Library","Adjusted P-value"]].to_dict("records") if len(dn_df) else []),
    }
    with open(os.path.join(OUT, "step2_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=str)
    print("[step2] 完成。")
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=str)[:1500])

if __name__ == "__main__":
    main()
