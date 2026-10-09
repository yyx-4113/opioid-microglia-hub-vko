# -*- coding: utf-8 -*-
"""
Step 3 — WGCNA-lite（含真实 TOM，BLAS 矩阵乘实现）
在阿片特异性签名基因集上构建共表达网络，软阈值 → 邻接 → TOM →
层次聚类 → 模块 → 模块特征基因(ME) ↔ 性状(吗啡/NP/对照)相关 →
关键模块内的模块内连通度最高者 = 枢纽候选。
输出：03_results/step3_module_assignment.csv, step3_module_trait.csv,
      step3_hub_candidates.csv, step3_summary.json
注：STRING 不可达，故"交互"以共表达 TOM 网络替代（v2 §6 接受的替代方案）。
"""
import os, json, sys
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.stats import pearsonr
import numpy.linalg as la
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "00_pipeline"))
import config as C

OUT = C.RESULTS
os.makedirs(OUT, exist_ok=True)

def load_logmat():
    fpkm = C.VERIFIED_GEO["GSE117320"]["files"]["fpkm"]
    df = pd.read_csv(fpkm, sep="\t", compression="gzip")
    df = df.set_index(df.columns[0])
    keep = [c for c in df.columns if True]
    # 仅保留带分组标记的样本列
    def grp(c):
        s = str(c).lower()
        if s.startswith("ctrl"): return "Control"
        if str(c).startswith("Mor"): return "Morphine"
        if any(k in s for k in ("snt","cci","sni")): return "NP"
        return None
    cols = [(c, grp(c)) for c in df.columns if grp(c) is not None]
    df = df[[c for c,_ in cols]].apply(pd.to_numeric, errors="coerce").dropna(how="all")
    df = df.loc[df.index.notna()]
    df = df[~df.index.astype(str).str.startswith("__")]
    if df.index.duplicated().any(): df = df.groupby(level=0).mean()
    log = np.log2(df.values.astype(float) + 1.0)
    groups = [g for _,g in cols]
    return df.index.values, log, np.array(groups), df.columns.values

def pick_power(corr):
    # scale-free topology fit: R^2 of log(p(k)) vs log(k) over powers
    # 仅 30 样本，网络不会完美 scale-free；在合理区间(4..14)取最小达到 fit>=0.8 的 β
    ks = []
    for b in range(4, 15):
        a = np.abs(corr)**b
        k = a.sum(1) - 1
        kk = k[k>0]; pk = np.bincount(kk.astype(int))
        deg = np.arange(len(pk)); p = pk/pk.sum()
        nz = p>0
        if nz.sum() < 3: continue
        x = np.log10(deg[nz]+1); y = np.log10(p[nz]+1e-12)
        if np.var(x) == 0: continue
        r = np.corrcoef(x, y)[0,1]
        ks.append((b, float(r*r)))
    if not ks: return 6, []
    fit_best = max(f[1] for f in ks)
    # 最小 β 达到 0.8，否则取 fit 最高者，否则默认 6
    cand = [b for b,f in ks if f >= 0.8]
    if cand: return min(cand), ks
    return (max(ks, key=lambda t: t[1])[0]), ks

def main():
    genes, log, groups, sample_cols = load_logmat()
    # 限制到 NP 激活签名（up+down，重构后主签名）做模块检测
    up = pd.read_csv(os.path.join(OUT,"step1_activation_up_full.csv"))["gene"].astype(str).tolist()
    dn = pd.read_csv(os.path.join(OUT,"step1_activation_down_full.csv"))["gene"].astype(str).tolist()
    sig = sorted(set(up)|set(dn))
    mask = np.isin(genes, sig)
    gsub = genes[mask]; L = log[mask]   # n x 30
    # 剔除零方差基因（跨样本恒定的 FPKM 会产生 NaN 相关，污染整个连通度求和）
    var = L.var(axis=1)
    keepv = np.isfinite(var) & (var > 1e-6)
    L = L[keepv]; gsub = gsub[keepv]
    print(f"[step3] 签名基因(去零方差后) {len(gsub)} / 全矩阵 {len(genes)}")
    # 基因-基因相关（NaN 兜底为 0）
    corr = np.corrcoef(L)
    corr = np.nan_to_num(corr, nan=0.0)
    np.fill_diagonal(corr, 1.0)
    corr = np.clip(corr, -0.999, 0.999)
    b, ks = pick_power(corr)
    fitstr = f"{max(x[1] for x in ks):.3f}" if ks else "n/a"
    print(f"[step3] 软阈值 β={b} (scale-free fit={fitstr})")
    adj = np.power(np.abs(corr), b)
    k = adj.sum(1) - 1.0
    # TOM（向量化）
    adj2 = adj @ adj
    num = adj + adj2
    denom = np.minimum.outer(k, k) + 1.0 - adj
    tom = num / np.maximum(denom, 1e-12)
    tom = np.clip(tom, 0.0, 1.0)   # TOM 边界钳制，避免 diss<0
    diss = 1.0 - tom
    np.fill_diagonal(diss, 0.0)
    # 层次聚类
    iu = np.triu_indices(len(gsub), 1)
    condensed = diss[iu]
    Z = linkage(condensed, method="average")
    # ---- Louvain 社区发现（稀疏共表达网络更稳健，避免层次切割碎裂成单例）----
    import networkx as nx
    from networkx.algorithms.community import louvain_communities
    gidx = {g: i for i, g in enumerate(gsub)}
    iu = np.triu_indices(len(gsub), 1)
    em = tom[iu] > 0.05
    G = nx.Graph()
    G.add_nodes_from(gsub.tolist())
    for src, tgt, w in zip(gsub[iu[0][em]], gsub[iu[1][em]], tom[iu[0][em], iu[1][em]]):
        G.add_edge(src, tgt, weight=float(w))
    comms = louvain_communities(G, weight="weight", resolution=1.0, seed=42)
    labels = np.zeros(len(gsub), int)
    for ci, c in enumerate(comms, 1):
        for g in c:
            labels[gidx[g]] = ci
    # 未分配（孤立点）给独立标签
    nxt = len(comms) + 1
    for i in range(len(gsub)):
        if labels[i] == 0:
            labels[i] = nxt; nxt += 1
    print(f"[step3] Louvain 模块数={len(set(labels))}  含边图节点={G.number_of_nodes()} 边={G.number_of_edges()}")
    # 模块特征基因（第一主成分，符号任意）
    expr = L.T  # 30 x n
    mods_uniq = sorted(set(labels.tolist()))
    me = {}
    for mod in mods_uniq:
        idx = labels == mod
        if idx.sum() < 2:
            me[mod] = expr[:, idx][:, 0] - expr[:, idx][:, 0].mean()
            continue
        sub = expr[:, idx]  # 30 x m
        u, s, vt = la.svd(sub - sub.mean(0), full_matrices=False)
        me[mod] = u[:, 0]
    # 性状向量（30 样本）
    trait = pd.DataFrame({
        "Morphine": (groups=="Morphine").astype(float),
        "NP": (groups=="NP").astype(float),
        "Control": (groups=="Control").astype(float),
    }, index=sample_cols)
    rows = []
    for mod in mods_uniq:
        pc1 = me[mod]
        rec = {"module": int(mod), "n_genes": int((labels==mod).sum())}
        for t in trait.columns:
            if np.std(pc1) == 0:
                rec[f"corr_{t}"] = 0.0; rec[f"p_{t}"] = 1.0; continue
            r, p = pearsonr(pc1, trait[t].values)
            rec[f"corr_{t}"] = round(float(r), 4); rec[f"p_{t}"] = float(p)
        rows.append(rec)
    mt = pd.DataFrame(rows).sort_values("corr_NP", key=lambda s: s.abs(), ascending=False)
    mt.to_csv(os.path.join(OUT, "step3_module_trait.csv"), index=False)
    # 关键模块：与 NP(神经损伤/小胶质激活臂) |corr|>=0.3 且模块规模>=5（避免单例噪声）
    # 注：吗啡臂为转录组静默(见 step1_morphine_null_report)，故关键模块以真实有信号的 NP 臂为准
    key_mods = mt[(mt["corr_NP"].abs() >= 0.3) & (mt["n_genes"] >= 5)]["module"].tolist()
    print(f"[step3] 关键模块(与NP |r|>=0.3, n>=5)={key_mods}")
    # 模块内连通度（排除自环）→ 枢纽候选
    hub_rows = []
    for mod in key_mods:
        idx = labels == mod
        sub_tom = tom[np.ix_(idx, idx)]
        diag = np.diag(sub_tom)
        intramod = sub_tom.sum(1) - diag
        sub_genes = gsub[idx]
        order = np.argsort(intramod)[::-1][:20]
        for rank, gi in enumerate(order):
            hub_rows.append({"module": int(mod), "gene": sub_genes[gi],
                             "intramodular_k": round(float(intramod[gi]), 4),
                             "rank_in_module": int(rank+1)})
    hubs = pd.DataFrame(hub_rows)
    hubs.to_csv(os.path.join(OUT, "step3_hub_candidates.csv"), index=False)
    pd.DataFrame({"gene": gsub, "module": labels}).to_csv(
        os.path.join(OUT, "step3_module_assignment.csv"), index=False)
    summary = {
        "beta": int(b), "scale_free_fit": round(float(max(x[1] for x in ks)), 3),
        "n_signature_genes": int(len(gsub)), "n_modules": int(len(np.unique(labels))),
        "key_modules": [int(m) for m in key_mods],
        "G3_input_candidates": int(len(hubs)),
        "top_hubs": hubs.head(20)[["module","gene","intramodular_k"]].to_dict("records") if len(hubs) else [],
    }
    json.dump(summary, open(os.path.join(OUT,"step3_summary.json"),"w"), indent=2, ensure_ascii=False)
    print("[step3] 完成。")
    print(json.dumps(summary, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
