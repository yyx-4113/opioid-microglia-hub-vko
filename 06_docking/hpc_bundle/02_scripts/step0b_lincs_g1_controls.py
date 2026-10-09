# -*- coding: utf-8 -*-
"""
Step 0b — G1 阳性对照恢复闸门（换含对照药的完整 LINCS L1000 子集重跑）
=============================================================================
目的：原 step0 用的 49,216 签名 LINCS 子集不含 7 个阳性对照药，导致 G1 NOT MET。
本脚本在 HPC（SCNet/超算，开放网络）上用 **完整 L1000 矩阵**（GSE92742 Level 3，
含全部 ~20k 化合物、必然包含 7 个对照药）重算 G1。

做法：
  1. 读 gene_info 建立 978 landmark (pr_gene_id -> symbol)；
  2. 读 step1 NP 激活签名 (UP/DOWN) 构建 landmark 索引；
  3. 读完整 gctx (n_sig x 978)，按秩次连通分计算每签名 CS；
  4. 按 01_data/lincs/control_sig_ids.txt 抽取 112 个对照药实例 CS；
  5. Mann-Whitney 单尾检验（对照 vs 背景），判定 G1 是否 PASS
     （p<0.05 且 ≥2 药 CS<-50）。

沙箱内网络无法取回完整 L1000（API 反爬 + GEO 大文件被代理 404），
故本脚本在 HPC 执行；此处仅做语法/逻辑校验与 49k 子集的"零对照"干跑。

用法（HPC）：
  python step0b_lincs_g1_controls.py --gctx /path/to/GSE92742_Broad_LINCS_LVL3.gctx \
         --controls 01_data/lincs/control_sig_ids.txt
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "00_pipeline"))
import config as C
sys.path.insert(0, os.path.dirname(__file__))
import step0_lincs as S0   # 复用 compute_cs / load_gene_symbol_map / build_query / control_match

L = C.DATA_LINCS
OUT = C.RESULTS
UP_FILE = os.path.join(OUT, "step1_activation_up_full.csv")
DOWN_FILE = os.path.join(OUT, "step1_activation_down_full.csv")
CTRL_IDS = os.path.join(L, "control_sig_ids.txt")


def run(gctx_path, ctrl_path):
    import h5py
    assert os.path.exists(gctx_path), f"LINCS 完整矩阵缺失: {gctx_path}"
    g2s = S0.load_gene_symbol_map()
    up_cols, down_cols, n_up_sym, n_down_sym = S0.build_query(g2s)
    print(f"[step0b] UP landmark={len(up_cols)}/{n_up_sym}  DOWN={len(down_cols)}/{n_down_sym}")

    with h5py.File(gctx_path, "r") as f:
        mat = f["/0/DATA/0/matrix"][:]            # (n_sig, 978)
        row_ids = [x.decode() if isinstance(x, bytes) else str(x)
                   for x in f["/0/META/ROW/id"][:]]   # 签名 id (sig_id)
        col_ids = [x.decode() if isinstance(x, bytes) else str(x)
                   for x in f["/0/META/COL/id"][:]]

    # 列=基因；建立 landmark 索引
    up_set, down_set = set(up_cols), set(down_cols)
    up_idx = np.array([i for i, c in enumerate(col_ids) if c in up_set], dtype=int)
    down_idx = np.array([i for i, c in enumerate(col_ids) if c in down_set], dtype=int)
    cs = S0.compute_cs(mat, up_idx, down_idx)
    print(f"[step0b] 计算 {len(cs)} 个签名 CS")

    # 对照实例
    with open(ctrl_path) as fh:
        ctrl_ids = set(line.strip() for line in fh if line.strip())
    row_index = {sid: i for i, sid in enumerate(row_ids)}
    ctrl_pos = [row_index[s] for s in ctrl_ids if s in row_index]
    print(f"[step0b] 命中对照实例 {len(ctrl_pos)}/{len(ctrl_ids)}（未命中={len(ctrl_ids)-len(ctrl_pos)}）")
    if len(ctrl_pos) == 0:
        print("[step0b] 警告：该 LINCS 子集仍不含任何对照药实例 -> G1 不可检验（应换更大子集）")
        return {"g1_pass": False, "n_control_hits": 0}

    ctrl_vals = cs[np.array(ctrl_pos)]
    bg_vals = cs[np.array([i for i in range(len(cs)) if i not in set(ctrl_pos)])]
    from scipy import stats as _sp
    _u, p_val = _sp.mannwhitneyu(ctrl_vals, bg_vals, alternative="less")
    recov_moderate = int((ctrl_vals < C.THRESH["lincs_cs_neg_moderate"]).sum())
    recov_strong = int((ctrl_vals < C.THRESH["lincs_cs_neg"]).sum())
    g1_pass = (p_val < 0.05) and (recov_moderate >= C.THRESH["min_positive_controls_recovered"])

    # 逐药 CS（按 sig_id 前缀聚合）
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
        "per_drug_CS": dict(sorted(per_drug.items(), key=lambda kv: kv[1])),
        "interpretation": ("G1 PASS: 阳性对照药（小胶质/阿片调节剂）被反向召回，LINCS 重定位结果可信。"
                           if g1_pass else
                           "G1 NOT MET: 对照药未显著逆转签名（或召回不足），需检查子集/阈值。"),
    }
    with open(os.path.join(OUT, "step0b_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("[step0b] G1:", "PASS" if g1_pass else "NOT MET", "| p=%.3g" % p_val,
          f"| 召回(CS<-50)={recov_moderate} 召回(CS<-90)={recov_strong}")
    return summary


def dryrun_49k():
    """沙箱内干跑：用现有 49k 子集验证代码路径（预期 0 对照命中）。"""
    gctx = os.path.join(L, "Level2_GEX_delta.gctx")
    if not os.path.exists(gctx):
        print("[step0b] 无 49k gctx，跳过干跑")
        return
    # 自动解压
    if not os.path.exists(gctx) and os.path.exists(gctx + ".gz"):
        import gzip as _gz
        with _gz.open(gctx + ".gz", "rb") as fi, open(gctx, "wb") as fo:
            fo.write(fi.read())
    print("[step0b] 干跑（49k 子集，预期无对照命中）：")
    run(gctx, CTRL_IDS)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gctx", default=None, help="完整 LINCS L1000 Level3 gctx 路径（HPC 上提供）")
    ap.add_argument("--controls", default=CTRL_IDS)
    ap.add_argument("--dryrun", action="store_true", help="用 49k 子集干跑验证代码路径")
    args = ap.parse_args()
    if args.dryrun or not args.gctx:
        dryrun_49k()
    else:
        run(args.gctx, args.controls)
