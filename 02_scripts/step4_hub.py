# -*- coding: utf-8 -*-
"""
Step 4 — PPI(共表达替代) + 机器学习三重收敛筛选枢纽基因
STRING 不可达，故用 Step3 的 TOM 共表达网络模块(模块内连通度)作"互作"代理；
ML 三重过滤(LASSO / RandomForest-Boruta式重要性 / SVM-RFE)在 M vs C 标签下
对签名基因排序，取 ≥2 法共同选中且落在关键模块枢纽集的基因 = 最终枢纽基因(3–6)。
输出：03_results/step4_* + step4_summary.json (含 G3 闸门)
"""
import os, json, sys
from collections import Counter
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "00_pipeline"))
import config as C

OUT = C.RESULTS
os.makedirs(OUT, exist_ok=True)

def load_logmat():
    fpkm = C.VERIFIED_GEO["GSE117320"]["files"]["fpkm"]
    df = pd.read_csv(fpkm, sep="\t", compression="gzip").set_index(df_cols0(fpkm))
    def grp(c):
        s = str(c).lower()
        if s.startswith("ctrl"): return "Control"
        if str(c).startswith("Mor"): return "Morphine"
        if any(k in s for k in ("snt","cci","sni")): return "NP"
        return None
    cols = [(c, grp(c)) for c in df.columns if grp(c) is not None]
    df = df[[c for c,_ in cols]].apply(pd.to_numeric, errors="coerce").dropna(how="all")
    df = df.loc[df.index.notna()]; df = df[~df.index.astype(str).str.startswith("__")]
    if df.index.duplicated().any(): df = df.groupby(level=0).mean()
    return df.index.values, np.log2(df.values.astype(float)+1.0), np.array([g for _,g in cols]), np.array([c for c,_ in cols])

def df_cols0(fpkm):
    return pd.read_csv(fpkm, sep="\t", compression="gzip", nrows=0).columns[0]

def main():
    genes, log, groups, sample_cols = load_logmat()
    # 限定特征池 = NP 激活签名基因（重构后主签名）
    up = pd.read_csv(os.path.join(OUT,"step1_activation_up_full.csv"))["gene"].astype(str).tolist()
    dn = pd.read_csv(os.path.join(OUT,"step1_activation_down_full.csv"))["gene"].astype(str).tolist()
    sig = set(up)|set(dn)
    mask = np.isin(genes, list(sig))
    gsub = genes[mask]
    # 仅取 NP(神经损伤) vs Control 样本（小胶质激活分类，重构后主签名对应臂）
    keep = np.isin(groups, ["NP","Control"])
    X = log[mask][:, keep].T   # (n_np+n_ctrl, n_sig)
    y = (groups[keep]=="NP").astype(int)
    print(f"[step4] 特征=签名基因 {X.shape[1]}  样本(NP/Control)={X.shape[0]} (y=1:{int(y.sum())})")

    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
    from sklearn.feature_selection import RFE
    Xs = StandardScaler().fit_transform(X)

    # 1) LASSO (真 L1：elasticnet + l1_ratio=1，sklearn>=1.8 原 penalty='l1' 已弃用且退化为 L2)
    lr = LogisticRegression(penalty="elasticnet", l1_ratio=1.0, C=0.5,
                            solver="saga", max_iter=8000)
    lr.fit(Xs, y)
    sel_lasso = set(gsub[lr.coef_[0] != 0])
    # 2) RandomForest 重要性 Top150
    rf = RandomForestClassifier(n_estimators=600, random_state=42, n_jobs=-1)
    rf.fit(Xs, y)
    imp = rf.feature_importances_
    top_rf = set(gsub[np.argsort(imp)[::-1][:150]])
    # 3) SVM-RFE
    svc = SVC(kernel="linear", C=1.0)
    rfe = RFE(svc, n_features_to_select=max(20, X.shape[1]//20), step=0.1)
    rfe.fit(Xs, y)
    sel_svm = set(gsub[rfe.support_])
    print(f"[step4] LASSO={len(sel_lasso)}  RFtop150={len(top_rf)}  SVM-RFE={len(sel_svm)}")

    # 收敛计数（≥2 法）
    cnt = Counter()
    for g in sel_lasso: cnt[g]+=1
    for g in top_rf: cnt[g]+=1
    for g in sel_svm: cnt[g]+=1
    conv = [g for g,c in cnt.items() if c>=2]
    print(f"[step4] 三重收敛(≥2法)基因数={len(conv)}")

    # 与 WGCNA 关键模块交集，并用"真实模块内连通度"排序（重算签名基因相关）
    hubs = pd.read_csv(os.path.join(OUT,"step3_hub_candidates.csv"))  # 参考
    wgcna_hubs = set(hubs["gene"])
    mod_assign = pd.read_csv(os.path.join(OUT,"step3_module_assignment.csv"))
    try:
        key_mods = json.load(open(os.path.join(OUT,"step3_summary.json")))["key_modules"]
    except Exception:
        key_mods = sorted(set(hubs["module"].tolist()))
    # 重算签名基因相关矩阵（n_sig x n_sig，与 gsub 同序）；零方差行→全 0（nan_to_num），其 intramodular_k=0 自然不入选
    corr2 = np.nan_to_num(np.corrcoef(log[mask]), nan=0.0)
    g2i = {g: i for i, g in enumerate(gsub)}
    mod_of = dict(zip(mod_assign["gene"], mod_assign["module"]))
    mod_members = {}
    for m in key_mods:
        mem = [g2i[g] for g in mod_assign[mod_assign["module"] == m]["gene"] if g in g2i]
        if len(mem) >= 3:
            mod_members[m] = mem
    def intramod(g):
        mi = g2i.get(g)
        if mi is None: return 0.0
        m = mod_of.get(g)
        if m not in mod_members: return 0.0
        others = [j for j in mod_members[m] if j != mi]
        return float(np.abs(corr2[mi, others]).sum())
    # 最终候选 = 收敛且落在关键模块
    final_cand = [g for g in conv if mod_of.get(g) in set(key_mods)]
    if len(final_cand) < 3:
        final_cand = list(conv)
    # 过滤明显 gene-model / 看家噪声
    noise = lambda g: any(g.startswith(p) for p in ("Gm","Rik","BC0","Btbd","4930","Gt("))
    final_cand = [g for g in final_cand if not noise(g)]
    final_cand.sort(key=intramod, reverse=True)
    final_genes = final_cand[:6]

    rec = pd.DataFrame([{
        "gene": g,
        "n_methods": cnt[g],
        "in_lasso": g in sel_lasso, "in_rf_top150": g in top_rf, "in_svm_rfe": g in sel_svm,
        "in_wgcna_hub": g in wgcna_hubs,
        "intramodular_k": round(intramod(g), 3),
        "in_activation_signature": g in sig,
    } for g in final_genes])
    rec.to_csv(os.path.join(OUT,"step4_final_hubs.csv"), index=False)

    summary = {
        "n_features": int(X.shape[1]), "n_lasso": int(len(sel_lasso)),
        "n_rf_top150": int(len(top_rf)), "n_svm_rfe": int(len(sel_svm)),
        "n_converged_ge2": int(len(conv)),
        "n_final_after_wgcna_intersect": int(len(final_cand)),
        "final_hub_genes": final_genes,
        "G3_pass": bool(len(final_genes) >= 3),
    }
    json.dump(summary, open(os.path.join(OUT,"step4_summary.json"),"w"), indent=2, ensure_ascii=False)
    print("[step4] 完成。")
    print(json.dumps(summary, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
