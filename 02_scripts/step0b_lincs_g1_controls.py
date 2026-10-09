# -*- coding: utf-8 -*-
"""
Step 0b — G1 阳性对照恢复闸门（用含对照药的完整 LINCS L1000 Level3 重跑）
=============================================================================
真实完整矩阵：GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz
（约 5GB 压缩；解压后为 1319138 签名 × 12328 基因 ≈ 65GB，故**不解压到磁盘**，
而是 gzip 解到内存 BytesIO，h5py 分块读矩阵算 CS，利用节点大内存规避磁盘瓶颈）。

维度修正（关键）：
  - 真实 Level3 INF 是 12328 列基因；L1000 实际只测定 978 landmark，
    gene_info.txt.gz 也仅含 978 landmark 的 pr_gene_id->symbol 映射。
  - 正确做法（与 Broad/cMAP 标准一致）：在 12328 列中**筛出 978 landmark 列**，
    仅对这些列算秩次连通分（N=978），避免推断的 11350 基因稀释 landmark 信号。
  - 矩阵方向自适应：gctx 标准布局为 (gene=ROW, sig=COL)，读取后按列(签名)分块并转置。

用法（节点）：
  python step0b_lincs_g1_controls.py --gctx /path/GSE92742_LVL3.gctx.gz \
         --controls 01_data/lincs/control_sig_ids.txt
  python step0b_lincs_g1_controls.py --dryrun   # 用现有 49k 子集验证代码路径
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "00_pipeline"))
import config as C
sys.path.insert(0, os.path.dirname(__file__))
import step0_lincs as S0

L = C.DATA_LINCS
OUT = C.RESULTS
UP_FILE = os.path.join(OUT, "step1_activation_up_full.csv")
DOWN_FILE = os.path.join(OUT, "step1_activation_down_full.csv")
CTRL_IDS = os.path.join(L, "control_sig_ids.txt")


def _decode(ids):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in ids]


def run(gctx_path, ctrl_path):
    import h5py, io, gzip
    assert os.path.exists(gctx_path), f"LINCS 矩阵缺失: {gctx_path}"
    # 支持直接读 .gz（解压到内存，避免 65GB 磁盘占用；节点 RAM 充足）
    if gctx_path.endswith(".gz"):
        print(f"[step0b] 解压 {gctx_path} 到内存 (gzip -> BytesIO)...")
        with gzip.open(gctx_path, "rb") as fh:
            data = fh.read()
        print(f"[step0b] 内存 HDF5 = {len(data)/1e9:.1f} GB")
        hf = h5py.File(io.BytesIO(data), "r")
    else:
        hf = h5py.File(gctx_path, "r")
    dset = hf["/0/DATA/0/matrix"]
    n_row, n_col = dset.shape
    row_ids = _decode(hf["/0/META/ROW/id"][:])
    col_ids = _decode(hf["/0/META/COL/id"][:])
    if len(row_ids) == n_row and len(col_ids) == n_col:
        gene_ids, sig_ids = row_ids, col_ids
        transpose = True     # dset = (gene, sig)
    elif len(row_ids) == n_col and len(col_ids) == n_row:
        gene_ids, sig_ids = col_ids, row_ids
        transpose = False    # dset = (sig, gene)
    else:
        raise SystemExit(f"[step0b] 无法判定 gctx 布局: shape={dset.shape}, "
                         f"row_ids={len(row_ids)}, col_ids={len(col_ids)}")
    n_sig = len(sig_ids)
    print(f"[step0b] 矩阵 {dset.shape} (gene×sig={transpose}); n_sig={n_sig}, n_gene={len(gene_ids)}")

    # 978 landmark 映射
    g2s = S0.load_gene_symbol_map()                    # gid -> symbol(upper)
    landmark_set = set(g2s.keys())
    gene_to_col = {str(g): i for i, g in enumerate(gene_ids)}
    land_cols = [gene_to_col[g] for g in landmark_set if g in gene_to_col]
    print(f"[step0b] 12328 列中命中 landmark 列 = {len(land_cols)}/978")
    assert len(land_cols) >= 50, "landmark 列命中过少，检查 gene_info 与 gctx 基因 id 一致性"
    land_col_arr = np.array(land_cols, dtype=int)

    # UP/DOWN 符号 -> landmark 列索引
    up = pd.read_csv(UP_FILE); down = pd.read_csv(DOWN_FILE)
    up_sym = set(str(s).upper() for s in up["gene"].dropna())
    down_sym = set(str(s).upper() for s in down["gene"].dropna())
    sym2col = {}
    for i, g in zip(land_cols, [gene_ids[c] for c in land_cols]):
        s = g2s.get(str(g))
        if s:
            sym2col[s] = i
    up_idx = np.array([sym2col[s] for s in up_sym if s in sym2col], dtype=int)
    down_idx = np.array([sym2col[s] for s in down_sym if s in sym2col], dtype=int)
    print(f"[step0b] 签名→landmark: UP={len(up_idx)}/{len(up_sym)}  DOWN={len(down_idx)}/{len(down_sym)}")
    assert len(up_idx) > 0 and len(down_idx) > 0, "UP/DOWN landmark 索引为空"

    # 分块读签名，筛 landmark 978 列，算 CS（避免一次性 65GB 加载）
    chunk = 20000
    cs_parts = []
    for st in range(0, n_sig, chunk):
        en = min(st + chunk, n_sig)
        if transpose:
            block = np.asarray(dset[:, st:en], dtype=np.float32).T   # (chunk, n_gene)
        else:
            block = np.asarray(dset[st:en, :], dtype=np.float32)      # (chunk, n_gene)
        block_lm = block[:, land_col_arr]                            # (chunk, 978)
        cs_parts.append(S0.compute_cs(block_lm, up_idx, down_idx))
        if (st // chunk) % 5 == 0:
            print(f"[step0b] 进度 {en}/{n_sig}")
    cs = np.concatenate(cs_parts)
    print(f"[step0b] 已计算 {len(cs)} 个签名 CS")

    # 对照实例
    with open(ctrl_path) as fh:
        ctrl_ids = set(line.strip() for line in fh if line.strip())
    row_index = {sid: i for i, sid in enumerate(sig_ids)}
    ctrl_pos = [row_index[s] for s in ctrl_ids if s in row_index]
    print(f"[step0b] 命中对照实例 {len(ctrl_pos)}/{len(ctrl_ids)}（未命中={len(ctrl_ids)-len(ctrl_pos)}）")
    if len(ctrl_pos) == 0:
        print("[step0b] 警告：该 LINCS 子集仍不含任何对照药实例 -> G1 不可检验")
        return {"g1_pass": False, "n_control_hits": 0}

    ctrl_vals = cs[np.array(ctrl_pos)]
    bg_vals = cs[np.array([i for i in range(len(cs)) if i not in set(ctrl_pos)])]
    from scipy import stats as _sp
    _u, p_val = _sp.mannwhitneyu(ctrl_vals, bg_vals, alternative="less")
    recov_moderate = int((ctrl_vals < C.THRESH["lincs_cs_neg_moderate"]).sum())
    recov_strong = int((ctrl_vals < C.THRESH["lincs_cs_neg"]).sum())
    g1_pass = (p_val < 0.05) and (recov_moderate >= C.THRESH["min_positive_controls_recovered"])

    per_drug = {}
    for s in ctrl_ids:
        if s in row_index:
            per_drug[s] = float(cs[row_index[s]])

    summary = {
        "g1_pass": bool(g1_pass),
        "n_control_hits": len(ctrl_pos),
        "n_control_total": len(ctrl_ids),
        "mwu_one_sided_p": float(p_val),
        "recovered_moderate_CS_lt_-50": recov_moderate,
        "recovered_strong_CS_lt_-90": recov_strong,
        "n_sig_total": int(n_sig),
        "n_landmark_cols": int(len(land_cols)),
        "up_landmark": int(len(up_idx)),
        "down_landmark": int(len(down_idx)),
        "per_drug_CS": dict(sorted(per_drug.items(), key=lambda kv: kv[1])),
        "interpretation": ("G1 PASS: 阳性对照药（小胶质/阿片调节剂）被反向召回，LINCS 重定位结果可信。"
                           if g1_pass else
                           "G1 NOT MET: 对照药未显著逆转签名（或召回不足），需检查子集/阈值。"),
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "step0b_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("[step0b] G1:", "PASS" if g1_pass else "NOT MET", "| p=%.3g" % p_val,
          f"| 召回(CS<-50)={recov_moderate} 召回(CS<-90)={recov_strong}")
    return summary


def dryrun_49k():
    """节点/沙箱干跑：用现有 49k 子集验证代码路径（预期 0 对照命中）。"""
    gctx = os.path.join(L, "Level2_GEX_delta.gctx")
    if not os.path.exists(gctx):
        print("[step0b] 无 49k gctx，跳过干跑")
        return
    print("[step0b] 干跑（49k 子集，预期无对照命中）：")
    run(gctx, CTRL_IDS)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gctx", default=None)
    ap.add_argument("--controls", default=CTRL_IDS)
    ap.add_argument("--dryrun", action="store_true")
    args = ap.parse_args()
    if args.dryrun or not args.gctx:
        dryrun_49k()
    else:
        run(args.gctx, args.controls)
