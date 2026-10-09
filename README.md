# opioid-microglia-hub-vko

**Reproducible computational package for:** *Opioid exposure is largely transcriptionally quiescent in spinal cord microglia: a multi-omics dissection of the neuropathic-pain activation hub universe and its RUNX1-anchored, repurposable-analgesia target space.*

Corresponding author: Yongxin Yang, B.M. — The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China. ORCID 0009-0004-9698-6552.

This repository preserves the analysis scripts, intermediate results, docking input/output package, and figures behind the PLOS ONE submission. No new primary (human/animal) data were generated; all inputs are public, de-identified omics repositories.

> **Tag for the submitted manuscript:** `v6.0.0` (matches the Data Availability Statement in the manuscript §6).

---

## 1. Data sources (public, not deposited here)

| Dataset | Accession | Use | Obtain from |
|---|---|---|---|
| Mouse spinal-cord microglia RNA-seq | GSE117320 (NCBI GEO) | Step 0–6 transcriptomics | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117320 |
| LINCS L1000 Phase-I | GSE92742 (NCBI GEO) | Step 0b/7/8 signature recovery & drug repurposing | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742 |
| Independent human microglia | GSE306403 (NCBI GEO) | Step 9 cross-cell-type replication | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306403 |

The raw downloaded inputs that were used to *produce* `03_results/` are **not deposited** in this repository (they are public and large). Re-runners should download them from the GEO links above and place them under `01_data/` (see `00_pipeline/` for expected layout). The deposited `03_results/`, `04_figures/`, and `06_docking/` are the authoritative products referenced throughout the manuscript.

## 2. File map

| Path | Contents |
|---|---|
| `00_pipeline/` | Pipeline orchestration and run scripts (`run_pipeline.sh`, etc.). |
| `02_scripts/` | All analysis scripts: `step0_lincs.py` … `step9_human_anchor.py`, `step0b_lincs_g1_controls.py`, `make_decoy_control.py`, `make_figures.py`, `make_fig9.py`, and `build_submission_v6.py` (regenerates the manuscript docx + cover letter + highlights + reporting checklist). |
| `03_results/` | Intermediate outputs of every step (CS tables, enrichment tables, hub lists, quality-gate summaries, `step6_summary.json`, `step8b_vina.csv`, `step8c_decoy_summary.json`). |
| `04_figures/` | Result figures (Fig 1–9). |
| `06_docking/` | Real AutoDock Vina docking input/receptor/ligand package + HPC execution scripts + `decoy_control/` (random-box & decoy-ligand specificity panels) and `06_docking/README.md`. |
| `07_submission/` | Manuscript (`manuscript_plosone_v6.docx`), cover letter, highlights, reporting checklist. |

## 3. Software versions

- **Python** 3.11+ (analysis, ligand prep, decoy control, figure generation).
- **AutoDock Vina 1.2.7** — the deposited, version-matched authoritative docking run (local sandbox; all `kcal/mol` in the manuscript come from this run, `06_docking/decoy_control/results_*/`).
- **AutoDock Vina 1.2.3** — HPC (SCNet) pilot run; headline ≈ −8.50 kcal/mol, superseded by the 1.2.7 authoritative set (headline JAK3–amlexanox = −8.53 kcal/mol).
- **Open Babel**, **Meeko**, **RDKit** — ligand protonation / PDBQT preparation.
- **R 4.x + DESeq2** — source transcriptomics (author count-based analysis of GSE117320; our `step1_deg.py` parses the authoritative DESeq2 output).
- **Enrichr ChEA_2016** (web API) — transcription-factor target enrichment.
- **PubChem PUG-REST** (web API) — compound identifier resolution.

## 4. Reproduction

```bash
# 1) (optional) obtain public inputs into 01_data/ per 00_pipeline layout
# 2) transcriptomics + hub + drug-repurposing pipeline
bash 00_pipeline/run_pipeline.sh          # runs step1..step9 in order
# 3) docking (needs AutoDock Vina in PATH)
bash 06_docking/run_vina.sh               # real-box 30 tasks
bash 06_docking/decoy_control/run_decoy_control.sh   # random-box + decoy specificity
# 4) regenerate the submission package (docx/md)
python 02_scripts/build_submission_v6.py
```

## 5. Known caveats (reported honestly in the manuscript)

- Morphine (opioid) arms are **largely transcriptionally quiescent** in microglia — no detectable DEG program (n = 7 morphine samples; DESeq2 is under-powered for subtle effects, stated as a limitation).
- Quality gate **G1 NOT MET** (LINCS positive-control recovery: 0 hits at CS < −90; weak directional signal only).
- RUNX1 association is **hypothesis-generating** (human LINCS space, CS = −0.28, below the −1.0 reversal-hit threshold; not a formal reversal hit; not direct mouse-microglia evidence).
- Docking specificity is **significant but modest** (random-box null p ≈ 0.019; decoy null p ≈ 0.020; 40 Å permissive box).
- A Zenodo archival snapshot with MANIFEST checksums will be issued on acceptance.

## 6. License

Code and documentation: MIT (see `LICENSE`). Result data tables: CC-BY-4.0 where reuse applies. Source omics: GEO terms of the original depositions.
