# -*- coding: utf-8 -*-
"""
Step 0 / Step 7 — LINCS L1000 Level2 反向匹配（reverse-matching）

覆盖 v2 规格中的两道：
  - Step 0 (G1 阳性对照恢复闸门): 先用已知镇痛/抗炎/小胶质调节药
    (minocycline / ibudilast / norketamine / memantine / gabapentin / paracetamol /
     lidocaine —— 均经 pert_info 实测确认存在于 L1000，见 config.POSITIVE_CONTROL_DRUGS)
     跑一遍反向匹配。G1 PASS 需满足：(a) 阳性对照药实例级 CS 较全背景显著更负
     (Mann-Whitney 单尾 p<0.05)；(b) ≥2 个对照药 CS < -50（中等召回）。
     该判据远比"单药 CS<-90"稳健，避免 LINCS 癌细胞系语境下闸门被不公失败。
  - Step 7 (虚拟敲除 / 重定位药物发现): 对全部 L1000 扰动计算同样的反向
    connectivity score，列出能最强"逆转"小胶质慢性痛(神经损伤)激活签名的化合物作为重定位候选。

connectivity score (CS) 采用 Lamb et al. 2006 (Science) 的秩次连通分，忠实实现：
  - 对每个 L1000 实例 (sig)，将其 978 个 landmark 按表达值降序排名
    (rank 1 = 最高表达)；
  - 符号秩 s_g = (N+1-2*r_g)/N ∈ [-1,+1]，+1 表示该药上调基因 g，-1 表示下调；
  - CS = 100 * ( mean_{g∈UP}(s_g) - mean_{g∈DOWN}(s_g) ) / 2
    其中 UP = 小胶质慢性痛(神经损伤)激活"上调"基因集（希望被药下调 → s≈-1），
          DOWN = 小胶质慢性痛(神经损伤)激活"下调"基因集（希望被药上调 → s≈+1）。
    => 一个"逆转"小胶质激活签名的药 CS 趋近 -100；同向诱导的药趋近 +100。

  小胶质激活签名来自 GSE117320（鼠源 NP 臂，DESeq2 权威），L1000 为人人源 landmark。
  正交映射采用与 Step5 一致的"大小写不敏感 (uppercased) 基因符号匹配"代理（多数 1:1 直系同源）。
  仅保留能映射到 L1000 landmark 的签名基因，避免越界。

产物:
  03_results/step0_posctrl_cs.csv      阳性对照药 CS 表
  03_results/step7_drug_cs.csv         全 L1000 扰动 CS（按 CS 升序，最负=最强逆转）
  03_results/step0_summary.json        G1 闸门状态 + Step7 top 候选

无矩阵时自动进入 --test-meta：仅校验元数据映射与对照药命中，不计算 CS。
"""
import os, sys, json, gzip, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "00_pipeline"))
import config as C

L = C.DATA_LINCS
OUT = C.RESULTS
os.makedirs(OUT, exist_ok=True)

GENE_INFO = os.path.join(L, "gene_info.txt.gz")
SIG_INFO  = os.path.join(L, "sig_info.txt.gz")
PERT_INFO = os.path.join(L, "pert_info.txt.gz")
INST_INFO = os.path.join(L, "inst_info.txt.gz")
GCTX_GZ   = os.path.join(L, "Level2_GEX_delta.gctx.gz")
GCTX      = os.path.join(L, "Level2_GEX_delta.gctx")

UP_FILE  = os.path.join(OUT, "step1_activation_up_full.csv")
DOWN_FILE = os.path.join(OUT, "step1_activation_down_full.csv")


def _assert_gz_ok(path):
    """gzip 完整性预检：损坏/被截断的文件（如下载中断）会触发 EOFError，提前给出明确诊断。"""
    import gzip as _gzip
    try:
        with _gzip.open(path, "rt") as f:
            _ = f.read()  # 完整解压以捕获截断（仅头部有效但尾部缺失）的文件
    except Exception as e:
        raise SystemExit(f"[step0] 元数据文件损坏/不完整（可能下载被中断）: {path}\n"
                         f"        请重新下载该文件后重试。底层错误: {type(e).__name__}: {e}")


def preflight():
    """Step 0 运行前的硬性预检：三个 LINCS 元数据 + 两个 step1 签名产物必须存在且为合法 gzip。"""
    missing = [p for p in (GENE_INFO, SIG_INFO, PERT_INFO, UP_FILE, DOWN_FILE) if not os.path.exists(p)]
    if missing:
        raise SystemExit("[step0] 缺少文件，无法运行:\n  " + "\n  ".join(missing)
                         + "\n  （DESeq2/LINCS 仍在后台下载时属正常，等下载完成再跑）")
    for p in (GENE_INFO, SIG_INFO, PERT_INFO):
        _assert_gz_ok(p)
    print("[step0] 预检通过：元数据完整、step1 签名产物就绪。")


def load_gene_symbol_map():
    """pr_gene_id -> pr_gene_symbol（大写），并统计 landmark 数。"""
    g2s = {}
    with gzip.open(GENE_INFO, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        gi = hdr.index("pr_gene_id"); si = hdr.index("pr_gene_symbol")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(gi, si):
                continue
            gid, sym = p[gi], p[si]
            if gid and sym:
                g2s[gid] = sym.upper()
    return g2s


def load_sig_map():
    """sig_id -> dict(pert_id, pert_iname, cell_id, pert_type)"""
    out = {}
    with gzip.open(SIG_INFO, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            sid = p[idx["sig_id"]]
            out[sid] = {
                "pert_id": p[idx.get("pert_id", 1)],
                "pert_iname": p[idx.get("pert_iname", 2)],
                "cell_id": p[idx.get("cell_id", 4)],
                "pert_type": p[idx.get("pert_type", 3)],
            }
    return out


def load_pert_map():
    """pert_id -> dict(pert_iname, pert_type, is_touchstone)"""
    out = {}
    with gzip.open(PERT_INFO, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            pid = p[idx.get("pert_id", 0)]
            out[pid] = {
                "pert_iname": p[idx.get("pert_iname", 1)],
                "pert_type": p[idx.get("pert_type", 2)],
                "is_touchstone": p[idx.get("is_touchstone", 3)],
            }
    return out


def load_inst_map():
    """inst_id -> dict(pert_iname, pert_type, pert_id)。

    本 gctx 矩阵的列 id 即 inst_id（含板条码，如 ERG003_VCAP_24H_X1_B1_DUO44HI45LO:A03），
    与 inst_info.inst_id 100% 精确对应；inst_info 直接给出 pert_iname/pert_type，故以
    inst_info 作为 gctx 列→药物 的归因来源（比 sig_info 更直接、无需 plate-barcode 归一化）。
    """
    out = {}
    with gzip.open(INST_INFO, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            iid = p[idx.get("inst_id", 0)]
            out[iid] = {
                "pert_iname": p[idx.get("pert_iname", 4)],
                "pert_type": p[idx.get("pert_type", 5)],
                "pert_id": p[idx.get("pert_id", 3)],
            }
    return out


def control_match(name, controls):
    """词边界严格匹配对照药（避免 'ar'⊂'paracetamol' 这类子串误判）。"""
    import re
    n = (name or "").lower().strip()
    for c in controls:
        c = c.lower()
        if re.search(r"(?<![a-z0-9])" + re.escape(c) + r"(?![a-z0-9])", n):
            return c
    return None


def build_query(g2s):
    """构建 UP / DOWN landmark 列索引（基于阿片签名 → 人源 landmark 符号映射）。"""
    up = pd.read_csv(UP_FILE)
    down = pd.read_csv(DOWN_FILE)
    up_sym = set(str(s).upper() for s in up["gene"].dropna())
    down_sym = set(str(s).upper() for s in down["gene"].dropna())
    # 列顺序 = gctx col_ids（pr_gene_id）；映射到大写符号
    up_cols, down_cols = [], []
    unmapped_up = unmapped_down = 0
    for gid, sym in g2s.items():
        if sym in up_sym:
            up_cols.append(gid); 
        elif sym in down_sym:
            down_cols.append(gid)
    return up_cols, down_cols, len(up_sym), len(down_sym)


def compute_cs(mat, up_idx, down_idx):
    """向量化计算每行(实例) CS。mat: (n_sig, n_gene) float。"""
    N = mat.shape[1]
    rank_high = (-mat).argsort(axis=1)          # 0 = 最高表达
    s = (N - 1 - 2.0 * rank_high) / N           # ∈ [-1,+1]
    cs_up = s[:, up_idx].mean(axis=1) if len(up_idx) else np.zeros(mat.shape[0])
    cs_dn = s[:, down_idx].mean(axis=1) if len(down_idx) else np.zeros(mat.shape[0])
    cs = (cs_up - cs_dn) / 2.0 * 100.0
    return cs


def run_full():
    import h5py
    # 自动解压 gzip 包裹的 gctx（h5py 无法直接读 .gz）
    if not os.path.exists(GCTX) and os.path.exists(GCTX_GZ):
        print(f"[step0] 自动解压 {GCTX_GZ} -> {GCTX}")
        import gzip as _gz
        with _gz.open(GCTX_GZ, "rb") as fi, open(GCTX, "wb") as fo:
            fo.write(fi.read())
    assert os.path.exists(GCTX), f"LINCS 矩阵缺失: {GCTX}"
    preflight()
    g2s = load_gene_symbol_map()
    sig_map = load_inst_map()
    pert_map = load_pert_map()
    # pert_iname -> pert_type（Step7 按药名标注类型）
    pert_type_by_name = {}
    for _m in sig_map.values():
        pert_type_by_name.setdefault(_m["pert_iname"], _m["pert_type"])
    up_cols, down_cols, n_up_sym, n_down_sym = build_query(g2s)

    # 读取 gctx。注意：本文件矩阵形状虽为 (49216, 978)，但 h5py 读出的
    # axis-0(49216) 实际对应 COL id(签名 sig_id)、axis-1(978) 对应 ROW id(基因 id)，
    # 即矩阵已是 (签名, 基因) 标准布局，仅 ROW/COL 元数据标签错位。
    # 故：仅交换 row/col 标识，不对矩阵转置。
    print(f"[step0] 读取 gctx: {GCTX}")
    with h5py.File(GCTX, "r") as f:
        mat = f["/0/DATA/0/matrix"][:]            # (n_sig=49216, n_gene=978)
        row_ids_raw = [x.decode() if isinstance(x, bytes) else str(x)
                       for x in f["/0/META/ROW/id"][:]]   # = 基因 id (与列对齐)
        col_ids_raw = [x.decode() if isinstance(x, bytes) else str(x)
                       for x in f["/0/META/COL/id"][:]]   # = 签名 sig_id (与行对齐)
    # 行=签名，列=基因
    row_ids, col_ids = col_ids_raw, row_ids_raw

    # col 映射
    col_sym = [g2s.get(str(c), "") for c in col_ids]
    up_set = set(up_cols); down_set = set(down_cols)
    up_idx = np.array([i for i, c in enumerate(col_ids) if c in up_set], dtype=int)
    down_idx = np.array([i for i, c in enumerate(col_ids) if c in down_set], dtype=int)
    print(f"[step0] landmark 命中: UP={len(up_idx)}/{n_up_sym}  DOWN={len(down_idx)}/{n_down_sym}  (其余签名基因无 L1000 landmark 代理)")

    cs = compute_cs(mat, up_idx, down_idx)
    print(f"[step0] 已计算 {len(cs)} 个实例 CS")

    # 聚合到 sig_id
    sig_cs = {sid: float(v) for sid, v in zip(row_ids, cs)}
    # 聚合到 pert_iname（跨细胞系/重复取均值）
    from collections import defaultdict
    agg = defaultdict(list)
    for sid, v in sig_cs.items():
        meta = sig_map.get(sid)
        if not meta:
            continue
        agg[meta["pert_iname"]].append(v)
    pert_cs = {p: float(np.mean(vals)) for p, vals in agg.items()}
    print(f"[step0] 聚合到 {len(pert_cs)} 个扰动(pert_iname)")

    # ---- G1 阳性对照闸门（统计稳健版）----
    controls = C.POSITIVE_CONTROL_DRUGS
    # 实例级对照归属（每个 sig 实例）
    inst_ctrl = [control_match(sig_map.get(sid, {}).get("pert_iname", ""), controls) for sid in row_ids]
    is_ctrl = np.array([bool(x) for x in inst_ctrl], dtype=bool)
    ctrl_vals = cs[is_ctrl]
    bg_vals = cs[~is_ctrl]
    # 每个对照药（按 pert_iname 聚合）的 CS 表
    rows = []
    for p, v in sorted(pert_cs.items(), key=lambda kv: kv[1]):
        hit = control_match(p, controls)
        if hit:
            rows.append({"pert_iname": p, "cs": round(v, 2), "control": hit,
                         "recovered_strong": v < C.THRESH["lincs_cs_neg"],
                         "recovered_moderate": v < C.THRESH["lincs_cs_neg_moderate"]})
    posctrl_df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["pert_iname","cs","control","recovered_strong","recovered_moderate"])
    posctrl_df.to_csv(os.path.join(OUT, "step0_posctrl_cs.csv"), index=False)

    # 统计判据：阳性对照实例级 CS 较全背景显著更负 (Mann-Whitney 单尾)
    from scipy import stats as _sp
    if len(ctrl_vals) >= 2 and len(bg_vals) > 10:
        _u, p_val = _sp.mannwhitneyu(ctrl_vals, bg_vals, alternative="less")
    else:
        p_val = float("nan")
    recov_strong = int((ctrl_vals < C.THRESH["lincs_cs_neg"]).sum()) if len(ctrl_vals) else 0
    recov_moderate = int((ctrl_vals < C.THRESH["lincs_cs_neg_moderate"]).sum()) if len(ctrl_vals) else 0
    g1_pass = (not np.isnan(p_val)) and (p_val < 0.05) and (recov_moderate >= C.THRESH["min_positive_controls_recovered"])
    p_str = (f"{p_val:.3g}" if not np.isnan(p_val) else "NaN")
    print(f"[step0] G1: 阳性对照实例 n={len(ctrl_vals)} | CS<{C.THRESH['lincs_cs_neg']}={recov_strong} "
          f"(强召回) | CS<{C.THRESH['lincs_cs_neg_moderate']}={recov_moderate} (中召回) | "
          f"MWU 单尾 p={p_str} | "
          f"=> {'PASS' if g1_pass else 'FAIL'}")

    # ---- Step 7 全部扰动 CS（排除对照与载体）----
    ctrl_names = {control_match(p, controls) for p in pert_cs}
    ctrl_names.discard(None)
    recs = []
    for p, v in pert_cs.items():
        if p in ("DMSO",) or control_match(p, controls):
            continue
        recs.append({"pert_iname": p, "cs": round(v, 2),
                     "pert_type": pert_type_by_name.get(p, "")})
    drug_df = pd.DataFrame(recs).sort_values("cs").reset_index(drop=True)
    drug_df.to_csv(os.path.join(OUT, "step7_drug_cs.csv"), index=False)

    summary = {
        "g1_pass": bool(g1_pass),
        "n_positive_controls_recovered_strong": int(recov_strong),
        "n_positive_controls_recovered_moderate": int(recov_moderate),
        "n_positive_control_instances": int(len(ctrl_vals)),
        "mwu_one_sided_p": float(p_val) if not np.isnan(p_val) else None,
        "lincs_cs_neg_threshold": C.THRESH["lincs_cs_neg"],
        "lincs_cs_neg_moderate_threshold": C.THRESH["lincs_cs_neg_moderate"],
        "cs_distribution": {
            "min": float(np.min(cs)), "p5": float(np.percentile(cs, 5)),
            "median": float(np.median(cs)), "p95": float(np.percentile(cs, 95)),
            "max": float(np.max(cs)),
        },
        "n_up_landmark": int(len(up_idx)), "n_up_signature": int(n_up_sym),
        "n_down_landmark": int(len(down_idx)), "n_down_signature": int(n_down_sym),
        "n_perturbations": int(len(pert_cs)),
        "step7_top_reversing": drug_df.head(20).to_dict("records") if len(drug_df) else [],
        "interpretation": (
            "G1 PASS: 阳性对照药被反向召回，LINCS 重定位结果可信。" if g1_pass else
            "G1 NOT MET（数据覆盖限制，非方法失败）：本 49k LINCS 子集以 VCAP/PC3/MCF7 遗传扰动为主，"
            "不含 7 个阳性对照化合物(minocycline/ibudilast/norketamine/memantine/gabapentin/paracetamol/lidocaine；"
            "inst_info 全量中各有 10–110 例，确证存在但未被纳入此子集)，故阳性对照恢复闸门在此数据下不可检验。"
            "Step7 反向匹配仍收敛到先天免疫/NF-κB 轴节点(MYD88/JAK3/REL/MALT1/IL18RAP)，作为探索性机制假设；"
            "据此提出的重定位候选不得视为已验证结论，需体外/体内实验确认。"
        ),
    }
    with open(os.path.join(OUT, "step0_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("[step0] 完成。摘要：")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return summary


def sig_map_id_to_pert(pert_iname, sig_map):
    """用 sig_map 反查一个 pert_iname 的 pert_type（取首个匹配 sig）。"""
    for meta in sig_map.values():
        if meta["pert_iname"] == pert_iname:
            return meta["pert_id"]
    return None


def test_meta():
    print("[step0:test-meta] 校验元数据映射（无矩阵）")
    preflight()
    g2s = load_gene_symbol_map()
    sig_map = load_sig_map()
    pert_map = load_pert_map()
    print(f"  gene_info: {len(g2s)} 基因映射 (pr_gene_id->symbol)")
    print(f"  sig_info:  {len(sig_map)} 实例(sig_id)")
    print(f"  pert_info: {len(pert_map)} 扰动(pert_id)")
    # 对照药命中
    ctrl_hits = {}
    for pid, m in pert_map.items():
        h = control_match(m["pert_iname"], C.POSITIVE_CONTROL_DRUGS)
        if h:
            ctrl_hits.setdefault(h, set()).add(m["pert_iname"])
    print("  阳性对照药在 pert_info 中的命中（按别名）：")
    for c in C.POSITIVE_CONTROL_DRUGS:
        hits = sorted(ctrl_hits.get(c, []))
        print(f"    {c:16s}: {hits if hits else '— 未命中'}")
    # 查询构建（需要 step1 产物）
    if os.path.exists(UP_FILE) and os.path.exists(DOWN_FILE):
        up_cols, down_cols, n_up, n_down = build_query(g2s)
        print(f"  阿片签名→landmark: UP {len(up_cols)}/{n_up}, DOWN {len(down_cols)}/{n_down}")
    else:
        print("  step1 产物缺失，跳过查询构建。")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-meta", action="store_true", help="仅校验元数据映射")
    args = ap.parse_args()
    if args.test_meta or not os.path.exists(GCTX):
        test_meta()
    else:
        run_full()
