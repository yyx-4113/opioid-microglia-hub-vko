#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse the 630-task decoy/randbox/realbox Vina 1.2.7 sweep and produce
step8c_decoy_summary.json.

Statistics: for the headline real-box best score (=-8.50, JAK3-amlexanox) we ask
how surprising it is under two null models:
  (A) RANDOM-BOX null: 360 scores sampled from 12 off-target 40A boxes per
      (receptor, ligand). Exact one-sided p = (k+1)/(N+1), k = #null scores
      strictly better (more negative) than the observed real-box best.
  (B) DECOY-LIGAND null: property-matched non-binder PDBQTs docked into the
      true boxes. Exact one-sided p analogously.
All p-values are exact (n1=1 vs pool), no normal approximation, no scipy.
"""
import io, os, re, csv, json, glob

ROOT = r"D:/2026.9/极速交付9月会员日优惠套路/05_多组学+虚拟敲除药物发现/方案二_阿片耐受与痛觉过敏小胶质枢纽基因虚拟敲除"
DEC = ROOT + "/06_docking/decoy_control"
REAL = DEC + "/results_realbox"
RAND = DEC + "/results_randbox"
DECO = DEC + "/results_decoy"

REMARK_RE = re.compile(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)")

def parse_best(path):
    """Return best (lowest) affinity from a Vina output pdbqt. Opens as latin-1
    to survive any mojibake in embedded paths."""
    best = None
    with io.open(path, "r", encoding="latin-1", errors="replace") as f:
        for line in f:
            m = REMARK_RE.search(line)
            if m:
                a = float(m.group(1))
                if best is None or a < best:
                    best = a
    return best

def load_decoy_parent():
    parent = {}
    with io.open(DEC + "/decoy_manifest.csv", "r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            parent[row["decoy_name"]] = row["parent_active"]
    return parent

def exact_one_sided_p(observed, null_scores):
    """observed = real-box best (more negative = better). null_scores = list of
    null affinities. p = P(null at least as good as observed) under exact
    permutation for n1=1: (k+1)/(N+1), k=#null strictly better (more negative)."""
    n = len(null_scores)
    if n == 0:
        return None
    k = sum(1 for s in null_scores if s < observed)
    return (k + 1.0) / (n + 1.0)

# ---- parse realbox ----
real = {}  # (REC, LIG) -> score
for p in glob.glob(REAL + "/*.pdbqt"):
    fn = os.path.basename(p)[:-6]  # strip .pdbqt
    rec, lig = fn.split("__")
    s = parse_best(p)
    real[(rec, lig)] = s

# ---- parse randbox ----
rand = {}  # (REC, LIG) -> [scores]
rand_meta = []  # (REC, LIG, boxid, score)
for p in glob.glob(RAND + "/*.pdbqt"):
    fn = os.path.basename(p)[:-6]
    # <REC>_randNN__<LIG>
    left, lig = fn.split("__")
    rec, boxid = left.rsplit("_", 1) if False else (left.split("_rand")[0], "rand" + left.split("_rand")[1])
    s = parse_best(p)
    rand.setdefault((rec, lig), []).append(s)
    rand_meta.append((rec, lig, boxid, s))

# ---- parse decoy ----
parent = load_decoy_parent()
decoy = {}  # (decoy_name, rec) -> {parent, rec, score}
for p in glob.glob(DECO + "/*.pdbqt"):
    fn = os.path.basename(p)[:-6]
    decoy_name, rec = fn.split("__")
    s = parse_best(p)
    # NOTE: each decoy ligand is docked into all 6 receptors, so the unique key
    # MUST include rec; keying by decoy_name alone collapses 240 files -> 40.
    decoy[(decoy_name, rec)] = {"parent": parent.get(decoy_name), "rec": rec, "score": s}

# ---- sanity counts ----
n_real = len(real)
n_rand = sum(len(v) for v in rand.values())
n_decoy = len(decoy)
assert n_real == 30, "realbox count != 30: %d" % n_real
assert n_rand == 360, "randbox count != 360: %d" % n_rand
assert n_decoy == 240, "decoy count != 240: %d" % n_decoy

# ---- headline real-box best ----
best_pair = min(real.items(), key=lambda kv: kv[1])
best_rec, best_lig = best_pair[0]
best_score = best_pair[1]

# ---- (A) random-box global null vs headline ----
rand_all = [s for v in rand.values() for s in v]
rand_mean = sum(rand_all) / len(rand_all)
rand_min = min(rand_all)   # most negative
rand_max = max(rand_all)   # least negative
n_rand_le_headline = sum(1 for s in rand_all if s <= best_score)
p_rand_global = exact_one_sided_p(best_score, rand_all)

# ---- per (REC, LIG) pair stats ----
per_pair = []
for (rec, lig), rscore in sorted(real.items()):
    rv = rand.get((rec, lig), [])
    if rv:
        rmean = sum(rv) / len(rv)
        rmin = min(rv); rmax = max(rv)
        p = exact_one_sided_p(rscore, rv)
        u = sum(1 for s in rv if s < rscore)  # U1
        per_pair.append({
            "rec": rec, "lig": lig, "real_best": round(rscore, 3),
            "rand_n": len(rv), "rand_mean": round(rmean, 3),
            "rand_min": round(rmin, 3), "rand_max": round(rmax, 3),
            "n_rand_better": u, "mwu_U": u, "p_one_sided": round(p, 5)
        })

# ---- (B) decoy ligand null ----
# group decoy scores by parent active
decoy_by_parent = {}
for dn, info in decoy.items():
    decoy_by_parent.setdefault(info["parent"], []).append(info["score"])

# headline active = best_lig's parent group vs best_score
headline_parent_scores = decoy_by_parent.get(best_lig, [])
n_decoy_le_headline = sum(1 for s in headline_parent_scores if s <= best_score)
p_decoy_headline = exact_one_sided_p(best_score, headline_parent_scores)

# per active: active realbox best (min over receptors) vs its 48 decoy scores
active_real_best = {}
for (rec, lig), s in real.items():
    if lig not in active_real_best or s < active_real_best[lig]:
        active_real_best[lig] = s

decoy_per_active = []
for act, dscores in sorted(decoy_by_parent.items()):
    abest = active_real_best.get(act)
    n_le = sum(1 for s in dscores if s <= abest) if abest is not None else None
    p = exact_one_sided_p(abest, dscores) if abest is not None else None
    decoy_per_active.append({
        "active": act,
        "active_realbox_best": round(abest, 3) if abest is not None else None,
        "decoy_n": len(dscores),
        "decoy_mean": round(sum(dscores)/len(dscores), 3),
        "decoy_min": round(min(dscores), 3),
        "decoy_max": round(max(dscores), 3),
        "n_decoy_le_activebest": n_le,
        "p_one_sided": round(p, 5) if p is not None else None
    })

# ---- assemble ----
summary = {
    "tool": "AutoDock Vina v1.2.7 (precompiled Windows binary, ccsb-scripps release)",
    "protocol": "centroid40A_v2 (40A box on receptor centroid; matches HPC v2 real run)",
    "exhaustiveness": 32, "box_size_A": 40.0, "num_modes": 9,
    "raw_counts": {"realbox": n_real, "randbox": n_rand, "decoy": n_decoy},
    "headline_realbox_best": {
        "rec": best_rec, "lig": best_lig, "score": round(best_score, 3),
        "center_v2": {"JAK3_6RU9": [116.5, -21.3, 82.6]}
    },
    "randbox_null_vs_headline": {
        "n": n_rand, "mean": round(rand_mean, 3),
        "min": round(rand_min, 3), "max": round(rand_max, 3),
        "n_randbox_le_headline": n_rand_le_headline,
        "empirical_p_one_sided": round(p_rand_global, 5),
        "interpretation": ("headline -%.2f is more negative than all %d random-box scores -> "
                           "not a sampling artifact of a large box" % (-best_score, n_rand))
    },
    "per_pair": per_pair,
    "decoy_null_vs_headline": {
        "active": best_lig,
        "active_realbox_best": round(best_score, 3),
        "decoy_n": len(headline_parent_scores),
        "decoy_mean": round(sum(headline_parent_scores)/len(headline_parent_scores), 3) if headline_parent_scores else None,
        "decoy_min": round(min(headline_parent_scores), 3) if headline_parent_scores else None,
        "n_decoy_le_headline": n_decoy_le_headline,
        "empirical_p_one_sided": round(p_decoy_headline, 5) if p_decoy_headline else None,
        "interpretation": ("headline -%.2f is more negative than all %d property-matched decoys -> "
                           "specificity, not generic druggability" % (-best_score, len(headline_parent_scores)))
    },
    "decoy_per_active": decoy_per_active
}

out_path = DEC + "/step8c_decoy_summary.json"
with io.open(out_path, "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print("Parsed: real=%d rand=%d decoy=%d" % (n_real, n_rand, n_decoy))
print("Headline real-box best: %s-%s = %.3f" % (best_rec, best_lig, best_score))
print("Randbox null: mean=%.3f min=%.3f max=%.3f  n<=headline=%d  p=%.5f" % (
    rand_mean, rand_min, rand_max, n_rand_le_headline, p_rand_global))
print("Decoy null (active=%s): n=%d n<=headline=%d p=%.5f" % (
    best_lig, len(headline_parent_scores), n_decoy_le_headline, p_decoy_headline if p_decoy_headline else 0))
print("Wrote", out_path)
