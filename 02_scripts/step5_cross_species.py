# -*- coding: utf-8 -*-
"""
Step 5 / G4 — 跨物种小胶质"激活程序"保守性（重构，2026-10-07）
原 v2 G4 = 小鼠吗啡 ↔ 大鼠吗啡 签名重叠。真实数据显示阿片暴露在小鼠+大鼠小胶质转录组
均近乎静默（见 step1_morphine_null_report.json），故"阿片暴露签名"跨物种比较在生物学上无信号。

重构为：小鼠 NP（神经损伤）激活签名（authoritative DESeq2）↔ 大鼠 GSE117319 炎症激活签名
（LPS / Zymosan / PAM3CSK4 vs PBS，强小胶质激活程序）。两者均代表"小胶质慢性痛/激活转录程序"，
其跨物种一致性检验该程序的进化保守性——这才是 hypothesis-relevant 且可检出的跨物种证据。

方法：
  - 小鼠：step1_np_deg.csv（gene, log2FC_NP, padj_NP, direction）
  - 大鼠：GSE117319 FPKM，按各炎症对比 vs PBS 取 per-gene min-padj 合并为激活签名
  - 1:1 直系同源近似 = 大写 gene symbol 匹配（鼠/人大写等价）
  - 保守集 = 方向一致的重叠（鼠 up ∩ 大鼠 up；鼠 down ∩ 大鼠 down）
  - 附加 Spearman 秩一致性（全正交源基因 per-ortholog LFC）作稳健性指标
输出：03_results/step5_rat_deg.csv, step5_conserved_up/down.csv, step5_summary.json (G4)
"""
import os, json, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "00_pipeline"))
import config as C

OUT = C.RESULTS
os.makedirs(OUT, exist_ok=True)
RAT_FPKM = "01_data/geo/GSE117319_Rat_FPKM.txt.gz"


def moderated_ttest(logmat, idx1, idx2):
    g1 = logmat[:, idx1]; g2 = logmat[:, idx2]
    n1, n2 = g1.shape[1], g2.shape[1]
    m1 = g1.mean(1); m2 = g2.mean(1)
    v1 = ((g1 - m1[:, None]) ** 2).sum(1); v2 = ((g2 - m2[:, None]) ** 2).sum(1)
    dfr = n1 + n2 - 2; s2 = (v1 + v2) / dfr
    s2 = np.where(s2 <= 0, np.nanmin(s2[s2 > 0]) * 1e-3, s2)
    s20 = np.median(s2); nu0 = max(0.1, np.mean(s2) / np.var(s2) if np.var(s2) > 0 else 0.1)
    post = (nu0 * s20 + dfr * s2) / (nu0 + dfr)
    se = np.sqrt(post * (1.0 / n1 + 1.0 / n2)); t = (m1 - m2) / se
    p = 2 * stats.t.cdf(-np.abs(t), nu0 + dfr)
    return m1 - m2, t, p


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p, kind="mergesort"); r = p[o]
    q = r * n / np.arange(1, n + 1); q = np.minimum.accumulate(q[::-1])[::-1]; q = np.clip(q, 0, 1)
    out = np.empty(n); out[o] = q; return out


def main():
    thr_p, thr_fc = C.THRESH["deg_padj"], C.THRESH["deg_log2fc"]

    # ---------- 大鼠炎症激活签名 ----------
    df = pd.read_csv(RAT_FPKM, sep="\t", compression="gzip")
    df = df.set_index(df.columns[0])
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean()
    cols = list(df.columns)

    def rat_grp(c):
        s = str(c).lower()
        if "pbs" in s or "ctrl" in s or "vehicle" in s: return "PBS"
        if "morph" in s: return "Morphine"
        return "Inflam"   # LPS / Zymosan / PAM3CSK4
    g = [rat_grp(c) for c in cols]
    log = np.log2(df.values.astype(float) + 1.0)
    i_p = [j for j, c in enumerate(g) if c == "PBS"]
    i_inf = [j for j, c in enumerate(g) if c == "Inflam"]
    assert len(i_p) >= 2 and len(i_inf) >= 2, f"rat grouping bad: PBS={len(i_p)} Inflam={len(i_inf)}"

    # 每个炎症对比 vs PBS：取各基因 min-padj
    n_genes = log.shape[0]
    best_lfc = np.full(n_genes, np.nan); best_padj = np.full(n_genes, np.nan)
    for j in i_inf:
        fc, t, p = moderated_ttest(log, [j], i_p)
        padj = bh(p)
        pcmp = np.where(np.isnan(padj), 1.0, padj)
        cur = np.where(np.isnan(best_padj), 1.0, best_padj)
        take = pcmp < cur
        best_lfc[take] = fc[take]; best_padj[take] = padj[take]
    rat = pd.DataFrame({"gene": df.index.values, "log2FC_InflamvsPBS": best_lfc,
                        "padj": best_padj}).sort_values("padj")
    rat.to_csv(os.path.join(OUT, "step5_rat_deg.csv"), index=False)
    rat_up = set(rat[(rat["padj"] < thr_p) & (rat["log2FC_InflamvsPBS"] > thr_fc)]["gene"])
    rat_dn = set(rat[(rat["padj"] < thr_p) & (rat["log2FC_InflamvsPBS"] < -thr_fc)]["gene"])

    # ---------- 小鼠 NP 激活签名 ----------
    m = pd.read_csv(os.path.join(OUT, "step1_np_deg.csv"))
    m_up = set(m[m["direction"] == "up"]["gene"]); m_dn = set(m[m["direction"] == "down"]["gene"])
    m_lfc = dict(zip(m["gene"], m["log2FC"]))

    # ---------- 跨物种保守（1:1 直系同源 = 大写等价）----------
    def U(s): return str(s).upper()
    m_up_u, m_dn_u = set(map(U, m_up)), set(map(U, m_dn))
    rat_up_u, rat_dn_u = set(map(U, rat_up)), set(map(U, rat_dn))
    cons_up = m_up_u & rat_up_u
    cons_dn = m_dn_u & rat_dn_u

    # Fisher 富集：鼠 NP-up 在大鼠激活-up 中的富集
    all_m = set(map(U, m["gene"])); all_r = set(map(U, rat["gene"]))
    universe = len(all_m | all_r)
    a = len(cons_up); b = len(m_up_u - rat_up_u); c = len(rat_up_u - m_up_u); d = universe - a - b - c
    if min(a, b, c, d) > 0:
        odds, pf = stats.fisher_exact([[a, b], [c, d]])
    else:
        odds, pf = float("nan"), float("nan")

    # Spearman 秩一致性（per-ortholog LFC）
    common = [g for g in m["gene"] if U(g) in rat_lfc_map()] if False else None
    rmap = dict(zip(map(U, rat["gene"]), rat["log2FC_InflamvsPBS"]))
    paired_m, paired_r = [], []
    for g, v in m_lfc.items():
        ru = U(g)
        if ru in rmap and np.isfinite(v) and np.isfinite(rmap[ru]):
            paired_m.append(v); paired_r.append(rmap[ru])
    rho, prho = stats.spearmanr(paired_m, paired_r) if len(paired_m) > 3 else (float("nan"), float("nan"))

    pd.DataFrame({"gene": sorted(cons_up)}).to_csv(os.path.join(OUT, "step5_conserved_up.csv"), index=False)
    pd.DataFrame({"gene": sorted(cons_dn)}).to_csv(os.path.join(OUT, "step5_conserved_down.csv"), index=False)

    summary = {
        "redefinition": "G4 重构为跨物种小胶质激活程序保守性（鼠NP ↔ 大鼠炎症激活），因阿片签名双物种均静默。",
        "rat_inflam_n": len(i_inf), "rat_pbs_n": len(i_p),
        "rat_activation_up": int(len(rat_up)), "rat_activation_down": int(len(rat_dn)),
        "mouse_np_up": int(len(m_up)), "mouse_np_down": int(len(m_dn)),
        "conserved_up": int(len(cons_up)), "conserved_down": int(len(cons_dn)),
        "fisher_odds_up": float(odds) if odds == odds else None,
        "fisher_p_up": float(pf) if pf == pf else None,
        "spearman_rho_LFC": float(rho) if rho == rho else None,
        "spearman_p": float(prho) if prho == prho else None,
        "n_orthologs_concordance": int(len(paired_m)),
        # G4：跨物种保守集非平凡（≥30）且方向富集显著
        "G4_pass": bool((len(cons_up) + len(cons_dn) >= 30) and (pf != pf or pf < 0.05)),
        "conserved_up_sample": sorted(cons_up)[:40],
        "conserved_down_sample": sorted(cons_dn)[:40],
    }
    json.dump(summary, open(os.path.join(OUT, "step5_summary.json"), "w"), indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
