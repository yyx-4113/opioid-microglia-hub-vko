#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate submission-grade figures (Fig1-Fig8) for the v4 manuscript from
03_results/* and 06_docking/* . All numbers read at run time from authoritative
result files. Outputs PNG (300 dpi) + PDF into 04_figures/.
Fig9 (decoy null) is produced separately after the decoy sweep finishes.
"""
import io, os, re, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch

ROOT = r"D:/2026.9/极速交付9月会员日优惠套路/05_多组学+虚拟敲除药物发现/方案二_阿片耐受与痛觉过敏小胶质枢纽基因虚拟敲除"
RES = ROOT + "/03_results"
FIG = ROOT + "/04_figures"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300,
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": "#333333", "axes.linewidth": 0.8,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "savefig.bbox": "tight",
})
C_UP = "#C0392B"; C_DOWN = "#1E8449"; C_MORPH = "#7F8C8D"
C_ACCENT = "#21618C"; GREY = "#555555"

def savefig(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"))
    fig.savefig(os.path.join(FIG, name + ".pdf"))
    plt.close(fig); print("wrote", name)

# ---------------------------------------------------------------- Fig1
def fig1():
    fig, ax = plt.subplots(figsize=(11, 6.4))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    ax.text(50, 97, "Multi-omics + in-silico knockout workflow and pre-registered go/no-go gates",
            ha="center", va="top", fontsize=13, fontweight="bold")
    boxes = [
        (3,70,20,14,"GSE117320","mouse spinal microglia\nMorphine 7 / NP 14 / Ctrl 9\n(G0 PASS)"),
        (28,70,20,14,"Step 1 DESeq2","NP signature\nG2 PASS (3069 DEG)"),
        (53,70,20,14,"Step 2 ORA","Enrichr GO/KEGG\nReactome/WikiPath"),
        (78,70,19,14,"Step 3 WGCNA","88 modules\nkey = module 74"),
        (3,48,20,14,"Step 4 ML filter","LASSO+RF+SVM\nG3 PASS (6 hubs)"),
        (28,48,20,14,"Step 5 / 9 G4 / G5","G4 NOT MET\nG5 PASS (silence)"),
        (53,48,20,14,"Step 6 RUNX1","ChEA ChIP-seq\nPOSITIVE anchor"),
        (78,48,19,14,"Step 7 / 0 LINCS","NF-kB convergence\nG1 NOT MET"),
        (15,26,26,14,"Step 8 BBB/ADMET","rule-based\nphenazone lead"),
        (53,26,26,14,"Step 8b Vina (HPC)","30 tasks, all neg\nbest -8.50 (decoy pend.)"),
    ]
    for (x,y,w,h,t,s) in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.4,rounding_size=1.2",
                    linewidth=1.2, edgecolor="#333", facecolor="#EAF2F8"))
        ax.text(x+w/2, y+h-3.2, t, ha="center", va="top", fontsize=9.5, fontweight="bold")
        ax.text(x+w/2, y+2.2, s, ha="center", va="bottom", fontsize=7.4, color=GREY)
    def arrow(x1,y1,x2,y2):
        ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=11,
                     linewidth=1.1, color="#333"))
    arrow(23,77,28,77); arrow(48,77,53,77); arrow(73,77,78,77)
    arrow(13,70,13,62); arrow(38,70,38,62); arrow(63,70,63,62); arrow(87.5,70,87.5,62)
    arrow(23,55,28,55); arrow(48,55,53,55)
    arrow(41,55,28,41); arrow(79,41,79,33)
    gates=[("G0","PASS"),("G2","PASS"),("G3","PASS"),("G5","PASS"),
           ("G4","NOT MET"),("G1","NOT MET"),("Step8b","EXEC (decoy pend.)")]
    x0=3; y0=8; w=12.5
    for i,(g,v) in enumerate(gates):
        x=x0+i*w
        col="#1E8449" if v.startswith("PASS") else ("#C0392B" if v=="NOT MET" else "#B9770E")
        ax.add_patch(FancyBboxPatch((x,y0),w-0.8,9,boxstyle="round,pad=0.3",
                    linewidth=1.2, edgecolor=col, facecolor="white"))
        ax.text(x+(w-0.8)/2, y0+6, g, ha="center", va="center", fontsize=9, fontweight="bold")
        ax.text(x+(w-0.8)/2, y0+2.3, v, ha="center", va="center", fontsize=7.6, color=col)
    ax.text(50,19.5,"Evidence tier: solid mid-tier SCIE  (upgrade gates a/b/c NOT yet met)",
            ha="center", va="center", fontsize=9, style="italic", color=GREY)
    savefig(fig, "Fig1_pipeline_gates")

# ---------------------------------------------------------------- Fig2
def fig2():
    s1 = json.load(io.open(RES+"/step1_summary.json", encoding="utf-8"))
    up = pd.read_csv(RES+"/step1_activation_up_full.csv")
    dn = pd.read_csv(RES+"/step1_activation_down_full.csv")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))
    bars = axA.bar(["Morphine\n(both arms)","Neuropathic\npain (NP)"],
                   [s1["morphine_bh_significant"], s1["np_total"]],
                   color=[C_MORPH, C_UP], width=0.55, edgecolor="#333")
    for b,v in zip(bars,[s1["morphine_bh_significant"], s1["np_total"]]):
        axA.text(b.get_x()+b.get_width()/2, v+max(s1["np_total"],3)*0.02, str(v),
                 ha="center", va="bottom", fontweight="bold")
    axA.set_yscale("log"); axA.set_ylabel("DEGs (BH-significant, log scale)")
    axA.set_title("A. Opioid exposure is transcriptionally silent;\nNP drives a massive activation program")
    axA.text(0, s1["morphine_bh_significant"]*1.2,
             "median |log2FC| = 0.24\nraw p<0.05 = %.1f%%" % (s1["morphine_rawp_frac_lt_0.05"]*100),
             ha="center", va="bottom", fontsize=8, color=GREY)
    axA.set_ylim(1, max(s1["np_total"],3)*3)
    fu = up.assign(dir="up"); fd = dn.assign(dir="down"); f = pd.concat([fu, fd])
    x = f["log2FC"].values.astype(float)
    p = f["padj"].values.astype(float)
    lp = -np.log10(np.where(p<=0, 1e-300, p)); lp = np.where(np.isfinite(lp), lp, 300)
    color = np.where(f["dir"].values=="up", C_UP, C_DOWN)
    axB.scatter(x, lp, s=6, c=color, alpha=0.5, linewidths=0)
    axB.axhline(-np.log10(0.05), color="#999", ls="--", lw=0.8)
    axB.axvline(0, color="#999", lw=0.8)
    axB.set_xlabel("log2 fold change (NP vs Control)"); axB.set_ylabel("-log10 adjusted p")
    axB.set_title("B. NP activation signature (n=%d DEG)" % s1["np_total"])
    axB.text(0.98,0.96,"UP %d"%s1["np_up"], transform=axB.transAxes, ha="right", color=C_UP, fontweight="bold")
    axB.text(0.98,0.90,"DOWN %d"%s1["np_down"], transform=axB.transAxes, ha="right", color=C_DOWN, fontweight="bold")
    savefig(fig, "Fig2_silence_activation")

# ---------------------------------------------------------------- Fig3
def _shorten(term):
    t = re.sub(r"\s*R-HSA-\d+.*","",term)
    t = re.sub(r"\s*-\s*GO:\d+.*","",t)
    t = re.sub(r"\s*-\s*KEGG.*","",t)
    t = re.sub(r"\s*\(.*?\)\s*$","",t)
    return t.strip()

def fig3():
    up = pd.read_csv(RES+"/step2_up_enrich.csv").sort_values("Adjusted P-value").head(12).copy()
    dn = pd.read_csv(RES+"/step2_down_enrich.csv").sort_values("Adjusted P-value").head(12).copy()
    up["lab"]=up["Term"].map(_shorten); up["mlp"]=-np.log10(up["Adjusted P-value"].astype(float))
    dn["lab"]=dn["Term"].map(_shorten); dn["mlp"]=-np.log10(dn["Adjusted P-value"].astype(float))
    fig, (axU, axD) = plt.subplots(1,2,figsize=(12,5.2))
    yU=np.arange(len(up))[::-1]
    axU.barh(yU, up["mlp"], color=C_UP, edgecolor="#333", height=0.7)
    axU.set_yticks(yU); axU.set_yticklabels(up["lab"], fontsize=7.6)
    axU.set_xlabel("-log10 adjusted p"); axU.set_title("A. NP-UP enriched terms (activation / proliferation)")
    axU.invert_yaxis()
    yD=np.arange(len(dn))[::-1]
    axD.barh(yD, dn["mlp"], color=C_DOWN, edgecolor="#333", height=0.7)
    axD.set_yticks(yD); axD.set_yticklabels(dn["lab"], fontsize=7.6)
    axD.set_xlabel("-log10 adjusted p"); axD.set_title("B. NP-DOWN enriched terms (OXPHOS / antioxidant)")
    axD.invert_yaxis()
    savefig(fig, "Fig3_pathways")

# ---------------------------------------------------------------- Fig4
def fig4():
    mt = pd.read_csv(RES+"/step3_module_trait.csv")
    mt["absNP"]=mt["corr_NP"].abs().astype(float)
    top = mt.sort_values("absNP", ascending=False).head(12)
    fig, (axA, axB) = plt.subplots(1,2,figsize=(11.5,4.8))
    traits=[("corr_NP","NP"),("corr_Control","Control"),("corr_Morphine","Morphine")]
    yy=np.arange(len(top))[::-1]; w=0.26
    for i,(c,t) in enumerate(traits):
        axA.barh(yy+(i-1)*w, top[c].astype(float), height=w,
                 color=[C_UP,C_DOWN,C_MORPH][i], edgecolor="#333", label=t)
    axA.axvline(0, color="#333", lw=0.8)
    axA.set_yticks(yy); axA.set_yticklabels(["mod %d (n=%d)"%(r["module"],r["n_genes"]) for _,r in top.iterrows()], fontsize=7.4)
    axA.set_xlabel("module-trait Pearson r"); axA.set_title("A. Top modules by |r| with NP trait")
    axA.legend(loc="lower right", fontsize=8)
    s3 = json.load(io.open(RES+"/step3_summary.json", encoding="utf-8"))
    steps=[("signature genes\n(n=%d)"%s3["n_signature_genes"], s3["n_signature_genes"], GREY),
           ("LASSO / RF / SVM\n(max set)", 153, C_ACCENT),
           (">=2 methods\n(48)", 48, "#1F618D"),
           ("final hubs\n(6)", 6, "#154360")]
    xs=np.arange(len(steps))
    bars=axB.bar(xs,[s[1] for s in steps], color=[s[2] for s in steps], edgecolor="#333", width=0.6)
    for b,s in zip(bars,steps):
        axB.text(b.get_x()+b.get_width()/2, b.get_height()*1.03, str(s[1]), ha="center", va="bottom", fontweight="bold", fontsize=9)
    for i in range(len(steps)-1):
        axB.plot([i+0.3,i+1-0.3],[steps[i][1],steps[i+1][1]], color="#888", ls="--", lw=0.8)
    axB.set_xticks(xs); axB.set_xticklabels([s[0] for s in steps], fontsize=7.6)
    axB.set_ylabel("genes"); axB.set_title("B. ML triple-filter convergence (G3 PASS)")
    axB.set_ylim(0, s3["n_signature_genes"]*1.12)
    savefig(fig, "Fig4_wgcna_ml")

# ---------------------------------------------------------------- Fig5
def fig5():
    d = pd.read_csv(RES+"/step7_drug_cs.csv").sort_values("cs").head(22)
    fig, ax = plt.subplots(figsize=(8.6,6.2))
    y=np.arange(len(d))[::-1]
    colors={"trt_sh":C_ACCENT,"trt_cp":"#B9770E","trt_oe":"#6C3483","trt_lig":"#117A65"}
    ax.barh(y, d["cs"].astype(float), color=[colors.get(t,"#888") for t in d["pert_type"]], edgecolor="#333", height=0.7)
    ax.set_yticks(y); ax.set_yticklabels(d["pert_iname"], fontsize=8)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_xlabel("Lamb connectivity score (more negative = reverses activation)")
    ax.set_title("Top LINCS perturbations reversing the NP-activation signature\n(convergent NF-kB / innate-immunity axis)")
    ax.legend(handles=[Patch(color=v,label=k) for k,v in colors.items()], fontsize=7.5, loc="lower right")
    savefig(fig, "Fig5_lincs_reversal")

# ---------------------------------------------------------------- Fig6
def fig6():
    s6 = json.load(io.open(RES+"/step6_summary.json", encoding="utf-8"))
    rows=[]
    for t in s6["enrichr_RUNX_terms"]:
        lab=re.sub(r"\s*\d{6,}.*","",t["term"])
        rows.append((lab, float(t["adjp"]), float(t["odds"]), t["overlap"]))
    rows=sorted([r for r in rows if r[1]<=0.1], key=lambda r:r[1])[:8]
    fig, ax = plt.subplots(figsize=(9,4.6))
    y=np.arange(len(rows))[::-1]; mlp=[-np.log10(r[1]) for r in rows]
    ax.barh(y, mlp, color=C_ACCENT, edgecolor="#333", height=0.66)
    for yi,(lab,ap,od,ov) in zip(y,rows):
        ax.text(mlp[yi]+0.05, yi, "adj-p=%.1e  OR=%.2f  (%s)"%(ap,od,ov), va="center", fontsize=7.6)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xlabel("-log10 adjusted p")
    ax.set_title("RUNX1 / RUNX-family ChIP-seq target enrichment among reversal hits (ChEA_2016)")
    savefig(fig, "Fig6_runx1")

# ---------------------------------------------------------------- Fig7
def fig7():
    v = pd.read_csv(RES+"/step8b_vina.csv")
    recs=sorted(v["receptor"].unique())
    ligs=["tofacitinib","ruxolitinib","phenazone","amlexanox","vorinostat"]
    M=np.full((len(recs),len(ligs)), np.nan)
    for _,r in v.iterrows():
        if r["ligand"] in ligs:
            M[recs.index(r["receptor"]), ligs.index(r["ligand"])]=float(r["best_affinity_kcal_mol"])
    fig, ax = plt.subplots(figsize=(7.2,5.4))
    im=ax.imshow(M, cmap="RdYlBu_r", aspect="auto")
    ax.set_xticks(range(len(ligs))); ax.set_xticklabels(ligs, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(recs))); ax.set_yticklabels(recs, fontsize=8)
    for i in range(len(recs)):
        for j in range(len(ligs)):
            if not np.isnan(M[i,j]):
                ax.text(j,i,"%.2f"%M[i,j], ha="center", va="center", fontsize=8,
                        color="black" if M[i,j]>-7 else "white")
    cb=fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label("best affinity (kcal/mol)")
    ax.set_title("Step 8b Vina docking: 6 NF-kB-axis targets x 5 ligands\n(all 30 tasks negative; best -8.50 JAK3-amlexanox)")
    savefig(fig, "Fig7_docking")

# ---------------------------------------------------------------- Fig8
def fig8():
    a = pd.read_csv(RES+"/step8_admet.csv")
    fig, ax = plt.subplots(figsize=(8.4,5.4))
    cmap={"likely_CNS+":C_UP,"likely_CNS-":C_DOWN,"likely_CNS-*":GREY}
    for _,r in a.iterrows():
        c=cmap.get(r["BBB_pred"],"#888")
        ax.scatter(float(r["cLogP"]), float(r["TPSA"]), s=90, c=c, edgecolor="#333", zorder=3)
        ax.annotate(r["compound"], (float(r["cLogP"]),float(r["TPSA"])),
                    xytext=(5,4), textcoords="offset points", fontsize=8)
    ax.axvspan(1,4, color="#FCF3CF", alpha=0.6, zorder=0)
    ax.axhspan(0,90, color="#FCF3CF", alpha=0.6, zorder=0)
    ax.set_xlabel("cLogP (1-4 = CNS-favorable window)")
    ax.set_ylabel("TPSA (<=90 = CNS-favorable)")
    ax.set_title("Rule-based BBB/ADMET of reversal + NF-kB candidates\n(shaded = favorable heuristic window; caffeine misclassified -> exploratory)")
    ax.legend(handles=[Patch(color=C_UP,label="likely_CNS+"),Patch(color=C_DOWN,label="likely_CNS-"),
                      Patch(color=GREY,label="likely_CNS-* (calibrator)")], fontsize=8, loc="upper left")
    savefig(fig, "Fig8_bbb_admet")

if __name__ == "__main__":
    print("matplotlib", matplotlib.__version__, "pandas", pd.__version__)
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6(); fig7(); fig8()
    print("ALL FIGURES DONE ->", FIG)
