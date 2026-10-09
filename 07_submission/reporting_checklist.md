# Reporting checklist & transparency statement (PLOS ONE submission)

**Study type:** Computational / in-silico re-analysis of public transcriptomics (no prospectively enrolled cohort).

**Applicable EQUATOR / discipline checklists.**
- Transcriptomics: MINSEQE framework (MINimal information about a high-throughput SEQuencing Experiment), extended to RNA-seq — satisfied via public repository citation (GEO), accession numbers (GSE117320, GSE92742), and deposited analysis scripts specifying all processing.
- No CONSORT (no RCT), no STROBE (no observational cohort), no PRISMA (no systematic review) applies.

**Methodological transparency (soundness-based review).**
- Pre-registered quality gates G0-G6 documented in manuscript §5 and `06_docking/README.md`.
- Software versions: DESeq2 (count-based calls from source authors), Enrichr ChEA_2016, AutoDock Vina 1.2.7 (local) / 1.2.3 (HPC), OpenBabel 3.2.1, Meeko 0.8.0, RDKit.
- Randomization / blinding: not applicable (secondary analysis of published data).
- Statistical tests: DESeq2 Wald test; Fisher exact for overlap; Mann-Whitney U for LINCS recovery; exact one-sided p for decoy panel — all justified in §2/§3.

[Upload this file as Supporting Information S1 Text if the submission system requests a checklist.]
