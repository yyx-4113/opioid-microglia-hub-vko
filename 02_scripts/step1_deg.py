# -*- coding: utf-8 -*-
"""
Step 1 — GSE117320 真实 DESeq2 解析（作者原始 count-based 分析，权威）

重要重构说明（2026-10-07，真实数据落地后）：
  原 v2 spec 假设"吗啡(阿片)暴露在小胶质产生稳健 DEG 签名"，真实 DESeq2 证伪该假设：
    - 吗啡双臂(Fixed/Escalating Morphine)经 BH 校正后仅 0–3 基因显著，合并原始 p<0.05 比例=4.0%≈随机期望 → 转录组层面**静默**。
    - 神经损伤(NP: SNT/CCI/SNI)臂产生 3069 个强 DEG（2085 up / 984 down），生物学合理且稳健。
  因此本步将**主签名改为 NP（神经损伤）小胶质激活签名**，作为"小胶质慢性痛激活枢纽宇宙"；
  吗啡近零结果单独诚实记录（step1_morphine_null_report.json），用于论文核心阴性发现。
  下游 step2–5 / step0 现消费 *_activation_* 文件（即 NP 签名），文件名已如实重命名。

文件结构（GSE117320_Mouse_DESeq2.xlsx, 单 sheet 'Sheet1'）：
  row0 = 条件标签（SNT 1dpi / CCI 2dpi / ... / Fixed Morphine / Escalating Morphine）
  row1 = 子表头（gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj）
  col0 = gene_id；其后每条件占 6 列：baseMean|log2FC|lfcSE|stat|pvalue|padj

输出：03_results/
  step1_np_deg.csv               全部 NP DEG（含统计）
  step1_activation_up_full.csv   上调用签名（全量，按 padj 排序）
  step1_activation_down_full.csv 下调用签名（全量）
  step1_activation_up.csv        截断 150（供 step0/7 LINCS）
  step1_activation_down.csv      截断 150
  step1_morphine_null_report.json 吗啡双物种近零证据
  step1_summary.json             摘要 + 闸门 G2（NP 签名≥20 → PASS）
"""
import os, sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "00_pipeline"))
import config as C

OUT = C.RESULTS
os.makedirs(OUT, exist_ok=True)


def parse_deseq2_blocks(path):
    """返回 (gene[], label->(lfc,padj,pvalue) 向量字典)。"""
    df = pd.read_excel(path, header=None)
    r0 = df.iloc[0]
    r1 = df.iloc[1]
    gene = df.iloc[2:, 0].astype(str).values
    blocks = []
    c = 1
    while c < df.shape[1]:
        lbl = r0[c]
        if pd.isna(lbl):
            c += 1
            continue
        subs = {r1[c + j]: j for j in range(6)}  # baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
        blocks.append((str(lbl), c, subs))
        c += 6
    vecs = {}
    for lbl, start, subs in blocks:
        lfc = pd.to_numeric(df.iloc[2:, start + subs["log2FoldChange"]], errors="coerce").values.astype(float)
        padj = pd.to_numeric(df.iloc[2:, start + subs["padj"]], errors="coerce").values.astype(float)
        pval = pd.to_numeric(df.iloc[2:, start + subs["pvalue"]], errors="coerce").values.astype(float)
        vecs[lbl] = (lfc, padj, pval)
    return gene, vecs


def arm_minpadj(labels, vecs):
    """对某臂内多个对比，逐基因取 padj 最小者对应的 lfc / pvalue。"""
    n = len(next(iter(vecs.values()))[0])
    best_lfc = np.full(n, np.nan)
    best_padj = np.full(n, np.nan)
    best_pv = np.full(n, np.nan)
    for L in labels:
        a, p, pv = vecs[L]
        pcmp = np.where(np.isnan(p), 1.0, p)
        cur = np.where(np.isnan(best_padj), 1.0, best_padj)
        take = pcmp < cur
        best_lfc[take] = a[take]
        best_padj[take] = p[take]
        best_pv[take] = pv[take]
    return best_lfc, best_padj, best_pv


def main():
    xlsx = C.VERIFIED_GEO["GSE117320"]["files"]["deseq2"]
    if not os.path.exists(xlsx) or os.path.getsize(xlsx) < 13_500_000:
        raise SystemExit("[step1] DESeq2 文件缺失或不完整（<%d B），中止。" % 13_500_000)

    gene, vecs = parse_deseq2_blocks(xlsx)
    all_labels = list(vecs.keys())
    morph_labels = [l for l in all_labels if "morph" in l.lower()]
    np_labels = [l for l in all_labels if any(k in l for k in ("SNT", "CCI", "SNI"))]
    print(f"[step1] 解析到 {len(all_labels)} 个对比：吗啡={morph_labels}，NP={np_labels}")

    lfc_m, padj_m, pv_m = arm_minpadj(morph_labels, vecs)
    lfc_n, padj_n, pv_n = arm_minpadj(np_labels, vecs)

    thr_p = C.THRESH["deg_padj"]
    thr_fc = C.THRESH["deg_log2fc"]

    # ---- 吗啡臂：诚实近零报告 ----
    m_sig = (padj_m < thr_p) & (np.abs(lfc_m) > thr_fc)
    m_rawp = (pv_m < thr_p) & (np.abs(lfc_m) > thr_fc)
    morph_report = {
        "labels": morph_labels,
        "n_genes": int(len(gene)),
        "bh_significant": int(m_sig.sum()),
        "bh_up": int(((padj_m < thr_p) & (lfc_m > thr_fc)).sum()),
        "bh_down": int(((padj_m < thr_p) & (lfc_m < -thr_fc)).sum()),
        "rawp_sig_fc": int(m_rawp.sum()),
        "rawp_any_lt_0.05_frac": float(np.nanmean(pv_m < 0.05)),
        "median_abs_lfc": float(np.nanmedian(np.abs(lfc_m))),
        "interpretation": "合并吗啡原始 p<0.05 比例≈5%（随机期望）→ 转录组层面静默；非 DEG 驱动。",
        "top_bh_genes": (
            [{"gene": gene[i], "lfc": float(lfc_m[i]), "padj": float(padj_m[i])}
             for i in np.argsort(padj_m)[:10] if np.isfinite(padj_m[i])]
        ),
    }
    with open(os.path.join(OUT, "step1_morphine_null_report.json"), "w") as f:
        json.dump(morph_report, f, indent=2)

    # ---- NP 臂：主签名（小胶质慢性痛激活宇宙）----
    n_up = (padj_n < thr_p) & (lfc_n > thr_fc)
    n_dn = (padj_n < thr_p) & (lfc_n < -thr_fc)
    np_df = pd.DataFrame({
        "gene": gene, "log2FC": lfc_n, "pvalue": pv_n, "padj": padj_n,
        "direction": np.select([n_up, n_dn], ["up", "down"], default="ns"),
    })
    np_df.to_csv(os.path.join(OUT, "step1_np_deg.csv"), index=False)

    up_full = np_df[n_up].sort_values("padj")
    dn_full = np_df[n_dn].sort_values("padj")
    up_full[["gene", "log2FC", "padj"]].to_csv(os.path.join(OUT, "step1_activation_up_full.csv"), index=False)
    dn_full[["gene", "log2FC", "padj"]].to_csv(os.path.join(OUT, "step1_activation_down_full.csv"), index=False)
    up_full[["gene", "log2FC", "padj"]].head(C.THRESH["signature_up"]).to_csv(
        os.path.join(OUT, "step1_activation_up.csv"), index=False)
    dn_full[["gene", "log2FC", "padj"]].head(C.THRESH["signature_down"]).to_csv(
        os.path.join(OUT, "step1_activation_down.csv"), index=False)

    summary = {
        "data_mode": "DESeq2(author,count-based,authoritative)",
        "n_genes": int(len(gene)),
        "np_up": int(n_up.sum()),
        "np_down": int(n_dn.sum()),
        "np_total": int((n_up | n_dn).sum()),
        "morphine_bh_significant": int(m_sig.sum()),
        "morphine_rawp_frac_lt_0.05": morph_report["rawp_any_lt_0.05_frac"],
        # G2 闸门：主签名（NP 激活）基因数 ≥ 20
        "G2_pass": bool((n_up | n_dn).sum() >= 20),
        "G2_note": "G2 现以 NP 激活签名为'小胶质慢性痛枢纽宇宙'（原吗啡签名假设已被真实数据证伪）。",
        "signature_up_n": int(min(len(up_full), C.THRESH["signature_up"])),
        "signature_down_n": int(min(len(dn_full), C.THRESH["signature_down"])),
        "top_np_up": up_full["gene"].head(20).tolist(),
        "top_np_down": dn_full["gene"].head(20).tolist(),
    }
    with open(os.path.join(OUT, "step1_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("[step1] 完成。摘要：")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
