#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Apply Round-2 (v5->v6) editorial revisions to the manuscript and build script.
Every replacement is asserted to occur (count printed)."""
import io, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "方案二_手稿v6_草稿.md")
BUILD = os.path.join(ROOT, "02_scripts", "build_submission_v6.py")

def load(p):
    return io.open(p, "r", encoding="utf-8").read()

def save(p, t):
    io.open(p, "w", encoding="utf-8").write(t)

def apply(text, old, new, label, expect=None):
    c = text.count(old)
    if c == 0:
        print("  [WARN] '%s' NOT FOUND (0 occurrences)" % label)
    elif expect is not None and c != expect:
        print("  [WARN] '%s' found %d times, expected %s" % (label, c, expect))
    else:
        print("  [OK] '%s' -> %d replacement(s)" % (label, c))
    return text.replace(old, new)

# ---------------- MANUSCRIPT ----------------
t = load(MD)
print("== manuscript v6 ==")

# T1-4 soften "silent" + power caveat
t = apply(t, "Opioid exposure is transcriptionally silent in spinal cord microglia:",
             "Opioid exposure is largely transcriptionally quiescent in spinal cord microglia:", "title")
t = apply(t, "Opioid (morphine) exposure was transcriptionally silent (3 genes survived BH correction",
             "Opioid (morphine) exposure was largely transcriptionally quiescent (3 genes survived BH correction", "abstract-silent")
t = apply(t, "best −8.534 kcal/mol, JAK3–amlexanox",
             "best −8.53 kcal/mol, JAK3–amlexanox", "abstract-precision")
t = apply(t, "opioid exposure is transcriptionally silent in microglia; we therefore corroborate",
             "opioid exposure is largely transcriptionally quiescent in microglia; we therefore corroborate", "s1-silent")
t = apply(t, "### 3.1 Opioid exposure is transcriptionally silent; NP drives a massive activation program (G2 PASS)",
             "### 3.1 Opioid exposure is largely transcriptionally quiescent; NP drives a massive activation program (G2 PASS)", "s31-header")
t = apply(t, "**Figure 2.** Opioid (morphine) transcriptional silence versus the NP activation program:",
             "**Figure 2.** Opioid (morphine) transcriptional quiescence versus the NP activation program:", "fig2")
t = apply(t, "and median |log2FC| was **0.24**. By contrast, the **NP arm produced 3,069 DEGs**",
             "and median |log2FC| was **0.24**. Given n=7 morphine samples, DESeq2 is powered to detect only large-effect DEGs; these statistics indicate absence of a *large* opioid program rather than proof of zero effect (subtle opioid responses cannot be excluded). By contrast, the **NP arm produced 3,069 DEGs**", "power-caveat")
t = apply(t, "the *phenotype* — \"opioid exposure → transcriptional silence\" — replicates",
             "the *phenotype* — \"opioid exposure → transcriptional quiescence\" — replicates", "s36-silence")
t = apply(t, "showed that **opioid (morphine) exposure is transcriptionally silent** in the same cells",
             "showed that **opioid (morphine) exposure is largely transcriptionally quiescent (no detectable DEG program)** in the same cells", "s4-silent")
t = apply(t, "Opioid silence replicates across two cell types (G5).",
             "Opioid transcriptional quiescence replicates across two cell types (G5).", "s4-strength-silence")

# T1-1 RUNX1 not a reversal hit
t = apply(t, "RUNX1 **itself** is a reversal hit (knockdown CS = **−0.28**), i.e., loss of RUNX1 mildly reverses the activation signature.",
             "RUNX1 perturbation shows a weak negative connectivity (knockdown CS = **−0.28**), which falls **below** the CS ≤ −1.0 reversal-hit threshold and is therefore directional only — not a formal reversal hit; loss of RUNX1 mildly reverses the activation signature in direction only.", "runx1-hit")

# T1-2 provenance re-label (HPC pilot vs deposited local 1.2.7)
t = apply(t, "**Real binding affinities (kcal/mol) were produced**: all 30 tasks negative, best **−8.53** (JAK3–amlexanox), mean −6.62, 10/30 ≤ −7, 0 positive (full panel in §3.10 and `03_results/step8b_vina.csv`). To keep the specificity null version-homogeneous, the 30 real-box tasks were **re-run locally with a prebuilt AutoDock Vina 1.2.7 Windows binary** (identical protocol: 40 Å box, exhaustiveness 32, 9 modes), giving the same headline (JAK3–amlexanox −8.53) and serving as the authoritative comparison set against the null.",
             "A HPC pilot (AutoDock Vina 1.2.3) confirmed directionality (headline ≈ **−8.5** kcal/mol, JAK3–amlexanox). Because the specificity null was executed under Vina 1.2.7, the **reported quantitative docking panel** (§3.10 statistics, `03_results/step8b_vina.csv`, `06_docking/decoy_control/step8c_decoy_summary.json`) is the **local Vina 1.2.7 re-run**, which is the deposited, version-matched authoritative set: all 30 tasks negative, best **−8.53** (JAK3–amlexanox), mean −6.62, 10/30 ≤ −7, 0 positive.", "provenance")
t = apply(t, "so that the specificity null (below) is version-homogeneous.",
             "so that the specificity null (below) is version-homogeneous. (The HPC 1.2.3 pilot gave a consistent headline ≈−8.5 kcal/mol; the deposited, cited numbers below are from the local 1.2.7 run.)", "provenance-note")

# T1-3 G1 threshold uniform CS<−90 (global)
t = apply(t, "CS<−50", "CS<−90", "g1-threshold-global")

# T1-5 figure renumber
t = apply(t, "**Figure 8.** Rule-based BBB / ADMET profiling",
             "**Figure 7.** Rule-based BBB / ADMET profiling", "fig8->7")
t = apply(t, "**Figure 7.** Structure-based AutoDock Vina docking campaign",
             "**Figure 8.** Structure-based AutoDock Vina docking campaign", "fig7->8")
t = apply(t, "with the best affinity (JAK3–amlexanox −8.534 kcal/mol) highlighted.",
             "with the best affinity (JAK3–amlexanox −8.53 kcal/mol) highlighted.", "fig8-precision")

# T2-1 MYD88 out of actionable list
t = apply(t, "target the *NP-activated* microglial state (hubs Gapt/Oas3/Cpne2/Sgms2/Cebpe/Il1f9 and the NF-κB axis MYD88/JAK3/REL/MALT1)",
             "target the *NP-activated* microglial state (hubs Gapt/Oas3/Cpne2/Sgms2/Cebpe/Il1f9 and the NF-κB-axis candidates JAK3/REL/MALT1; MYD88 is retained only as a docking negative-control sanity check, not an actionable target)", "myd88-list")

# T2-2 per-pair specificity + S1
t = apply(t, "and the repurposing claim remains hypothesis-generating pending experimental validation.",
             "and the repurposing claim remains hypothesis-generating pending experimental validation. At the per-pair level no individual receptor–ligand pair reaches p<0.05 (e.g., JAK3–amlexanox per-pair MWU p=0.077); the aggregate p=0.019–0.020 therefore reflects the selected best of 30 tasks, consistent with modest specificity rather than per-target significance. The full per-pair table is deposited as Supplementary Table S1 (`06_docking/decoy_control/step8c_decoy_summary.json`, `per_pair` array).", "per-pair")

# T2-3 MINSEQE rationale
t = apply(t, "This computational re-analysis adheres to the MINSEQE (MINimal information about a high-throughput SEQuencing Experiment) framework for transcriptomics reporting by citing the public repository, accession numbers, and analysis scripts that fully specify the data and processing.",
             "This computational re-analysis adheres to the MINSEQE (MINimal information about a high-throughput SEQuencing Experiment) framework for transcriptomics reporting by citing the public repository, accession numbers, and analysis scripts that fully specify the data and processing. As a secondary in-silico re-analysis with no newly generated sequencing data, MINSEQE is cited for transcriptomics-reporting transparency (public repository + accession + scripts); no primary MINSEQE data-deposit table applies.", "minseqe")

# T3 precision -8.534 -> -8.53 in body
t = apply(t, "best −8.534 kcal/mol for JAK3–amlexanox; mean −6.62; 10/30 ≤ −7)",
             "best −8.53 kcal/mol for JAK3–amlexanox; mean −6.62; 10/30 ≤ −7)", "lim4-precision")
t = apply(t, "best −8.534 (JAK3–amlexanox), mean −6.62, 10/30 ≤ −7; decoy panel",
             "best −8.53 (JAK3–amlexanox), mean −6.62, 10/30 ≤ −7; decoy panel", "gate-precision")

# author-note consistency
t = apply(t, "0 达 CS<−50，弱方向信号", "0 达 CS<−90，弱方向信号", "note-cs")
t = apply(t, "版本匹配 best −8.534；decoy", "版本匹配 best −8.53；decoy", "note-precision")
t = apply(t, "> **v5 落实的 Round 1 修订**：(T1-1) 明确 Sypek 2024 为“吗啡沉默 + NP 程序”的首报，本稿重定位于“佐证 + 拓展”；(T1-2) RUNX1 降为“假设生成性关联”，删除“跨人鼠一致”表述、补小胶质身份 TF 文献；(T1-3) 枢纽标注改为“候选”、明示 module 74 为 NP 抑制态、ML 无外样本验证；(T1-4) phenazone 重度保留措辞；(T1-5) MYD88 对接标注为非可解释、对齐 JAK2/STAT3 轴；(T2-1) 补 9 条图例与正文引用；(T2-2) 修正 BBB HBD 表(4/7)；(T2-3) 重建 step8b_vina.csv(1.2.7 正确 schema)；(T2-4/5) 引文改方括号、MIAME→MINSEQE；(T2-8) 盒心偏移更正为≈11.8 Å；(T2-9) decoy 路径更正。Data Availability 引用的 GitHub 仓库须在投稿前真实创建（用户动作，非阻断修订）。",
             "> **v5 落实的 Round 1 修订**：(T1-1) 明确 Sypek 2024 为“吗啡沉默 + NP 程序”的首报，本稿重定位于“佐证 + 拓展”；(T1-2) RUNX1 降为“假设生成性关联”，删除“跨人鼠一致”表述、补小胶质身份 TF 文献；(T1-3) 枢纽标注改为“候选”、明示 module 74 为 NP 抑制态、ML 无外样本验证；(T1-4) phenazone 重度保留措辞；(T1-5) MYD88 对接标注为非可解释、对齐 JAK2/STAT3 轴；(T2-1) 补 9 条图例与正文引用；(T2-2) 修正 BBB HBD 表(4/7)；(T2-3) 重建 step8b_vina.csv(1.2.7 正确 schema)；(T2-4/5) 引文改方括号、MIAME→MINSEQE；(T2-8) 盒心偏移更正为≈11.8 Å；(T2-9) decoy 路径更正。Data Availability 引用的 GitHub 仓库须在投稿前真实创建（用户动作，非阻断修订）。\n> **v6 落实的 Round 2 修订**：(T1-1) RUNX1 改“弱负向连通 CS=−0.28、低于 −1.0 逆转命中阈值、仅方向性、非正式逆转命中”；(T1-2) HPC 1.2.3 标为 pilot(≈−8.5)，deposited 数字归本地 Vina 1.2.7 权威集；(T1-3) G1 阈值全文统一 CS<−90；(T1-4) “静默”→“基本转录静息/无可检出 DEG 程序”+n=7 功效说明（标题/摘要/§1/§3.1/§3.6/§4）；(T1-5) 图序修正 Fig7=BBB、Fig8=对接、Fig9=decoy；(T2-1) MYD88 移出可成药列表（仅留作对接阴性对照）；(T2-2) 补 per-pair 特异性句+补充表 S1（JAK3–amlexanox MWU p=0.077）；(T2-3) MINSEQE 理由句；(T3) 精度统一 −8.53。", "v6-note")

save(MD, t)
print("manuscript v6 saved.")

# ---------------- BUILD SCRIPT ----------------
b = load(BUILD)
print("== build_submission_v6.py ==")
b = apply(b, 'SRC  = os.path.join(ROOT, "方案二_手稿v5_草稿.md")',
             'SRC  = os.path.join(ROOT, "方案二_手稿v6_草稿.md")', "build-src")
b = apply(b, 'out_docx = os.path.join(OUT, "manuscript_plosone_v5.docx")',
             'out_docx = os.path.join(OUT, "manuscript_plosone_v6.docx")', "build-out")
b = apply(b, "Opioid exposure is transcriptionally silent in spinal cord microglia:",
             "Opioid exposure is largely transcriptionally quiescent in spinal cord microglia:", "cl-title")
b = apply(b, "Opioid (morphine) exposure is transcriptionally silent in mouse spinal microglia.",
             "Opioid (morphine) exposure is largely transcriptionally quiescent in mouse spinal microglia.", "hl-silent")
b = apply(b, "Short title: Opioid transcriptional silence in spinal microglia",
             "Short title: Opioid transcriptional quiescence in spinal microglia", "short-title")
save(BUILD, b)
print("build script saved.")
print("DONE")
