#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_fig9.py — Fig9: decoy / random-box NULL distribution for the Step 8b
specificity control (AutoDock Vina v1.2.7, local, 40A centroid protocol).

Reads the raw Vina output pdbqt files (no hardcoded numbers) and the
step8c_decoy_summary.json produced by analyze_decoy_control.py, then draws a
two-panel histogram:
  (A) random-box null (360 off-target 40A centers) with headline real-box best
  (B) decoy-ligand null (240 property-matched non-binders) with the same line
Each panel annotates the exact one-sided empirical p-value.

Usage:
  python make_fig9.py
Outputs: 04_figures/Fig9_decoy_null.png  +  .pdf
"""
import io, os, re, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEC = os.path.join(ROOT, "06_docking", "decoy_control")
REAL = os.path.join(DEC, "results_realbox")
RAND = os.path.join(DEC, "results_randbox")
DECO = os.path.join(DEC, "results_decoy")
SUMMARY = os.path.join(DEC, "step8c_decoy_summary.json")
OUTDIR = os.path.join(ROOT, "04_figures")
os.makedirs(OUTDIR, exist_ok=True)

REMARK_RE = re.compile(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)")

def parse_best(path):
    best = None
    with io.open(path, "r", encoding="latin-1", errors="replace") as f:
        for line in f:
            m = REMARK_RE.search(line)
            if m:
                a = float(m.group(1))
                if best is None or a < best:
                    best = a
    return best

def collect(d):
    vals = []
    for p in glob.glob(os.path.join(d, "*.pdbqt")):
        s = parse_best(p)
        if s is not None:
            vals.append(s)
    return vals

real = collect(REAL)
rand = collect(RAND)
deco = collect(DECO)

real = np.array(real) if real else np.array([])
rand = np.array(rand) if rand else np.array([])
deco = np.array(deco) if deco else np.array([])

headline = float(real.min()) if real.size else None

# load summary p-values if present
p_rand = p_deco = None
if os.path.isfile(SUMMARY):
    with io.open(SUMMARY, "r", encoding="utf-8") as f:
        S = json.load(f)
    p_rand = S.get("randbox_null_vs_headline", {}).get("empirical_p_one_sided")
    p_deco = S.get("decoy_null_vs_headline", {}).get("empirical_p_one_sided")

headline_str = f"{headline:.3f}" if headline is not None else "NA"
print(f"[make_fig9] real={real.size} rand={rand.size} decoy={deco.size} "
      f"headline={headline_str} "
      f"p_rand={p_rand} p_decoy={p_deco}")

# ---- plot ----
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12,
                     "axes.labelsize": 11, "figure.dpi": 300})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

bins = np.linspace(-9.5, -4.5, 26)

# Panel A: random-box null
ax = axes[0]
if rand.size:
    ax.hist(rand, bins=bins, color="#4C72B0", alpha=0.85, edgecolor="white", linewidth=0.4)
if headline is not None:
    ax.axvline(headline, color="#C44E52", linewidth=2.2, linestyle="--",
               label=f"real-box best = {headline:.2f}")
    cnt = int((rand <= headline).sum())
    txt = f"random-box null: n={rand.size}\n#≤ headline = {cnt}"
    if p_rand is not None:
        txt += f"\nexact one-sided p = {p_rand:.4f}"
    ax.set_title("A. Random-box null (off-target 40Å centers)")
else:
    ax.set_title("A. Random-box null (pending)")
ax.set_xlabel("Best Vina affinity (kcal/mol, more negative = stronger)")
ax.set_ylabel("Count")
if rand.size and headline is not None:
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)

# Panel B: decoy-ligand null
ax = axes[1]
if deco.size:
    ax.hist(deco, bins=bins, color="#55A868", alpha=0.85, edgecolor="white", linewidth=0.4)
if headline is not None:
    ax.axvline(headline, color="#C44E52", linewidth=2.2, linestyle="--",
               label=f"real-box best = {headline:.2f}")
    cnt = int((deco <= headline).sum())
    txt = f"decoy-ligand null: n={deco.size}\n#≤ headline = {cnt}"
    if p_deco is not None:
        txt += f"\nexact one-sided p = {p_deco:.4f}"
    # place text box
    ax.text(0.97, 0.95, txt, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="#999", alpha=0.9))
    ax.set_title("B. Decoy-ligand null (property-matched non-binders)")
else:
    ax.set_title("B. Decoy-ligand null (pending)")
ax.set_xlabel("Best Vina affinity (kcal/mol, more negative = stronger)")
if deco.size and headline is not None:
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)

fig.suptitle("Fig9. Step 8b specificity control: docking affinities under null models",
             fontsize=12, y=1.01)
fig.tight_layout()
png = os.path.join(OUTDIR, "Fig9_decoy_null.png")
pdf = os.path.join(OUTDIR, "Fig9_decoy_null.pdf")
fig.savefig(png, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"[make_fig9] wrote {png} and {pdf}")
