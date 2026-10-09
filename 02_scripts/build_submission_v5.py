#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_submission.py — assemble 07_submission/ for the v4 manuscript.

Produces (Times New Roman 12pt):
  (1) manuscript_v4.docx          — converted from 方案二_手稿v4_草稿.md
                                    (author-note blockquotes excluded; tables rendered;
                                     figures referenced, PNGs copied to figures/)
  (2) cover_letter.docx + .md
  (3) highlights.docx + .md
  (4) response_letter_template.docx + .md

Usage: python build_submission.py
"""
import io, os, re, glob, shutil
from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, "方案二_手稿v5_草稿.md")
FIGS = os.path.join(ROOT, "04_figures")
OUT  = os.path.join(ROOT, "07_submission")
FIGOUT = os.path.join(OUT, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIGOUT, exist_ok=True)

FONT = "Times New Roman"

# ---------- helpers ----------
def set_base_style(doc):
    st = doc.styles["Normal"]
    st.font.name = FONT
    st.font.size = Pt(12)
    st.element.rPr.rFonts.set(__import__("docx").oxml.ns.qn("w:eastAsia"), FONT)

def add_heading(doc, text, level):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.name = FONT
        run.font.size = Pt({1: 15, 2: 13, 3: 12}.get(level, 12))
        run.font.color.rgb = RGBColor(0, 0, 0)
    return h

def add_para(doc, text, italic=False, bold=False, size=12, color=None, align=None):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.font.name = FONT
    r.font.size = Pt(size)
    r.italic = italic
    r.bold = bold
    if color is not None:
        r.font.color.rgb = color
    if align is not None:
        p.alignment = align
    return p

def inject_title_page(doc):
    """PLOS ONE requires author/affiliation/email/ORCID on the title page."""
    title = io.open(SRC, "r", encoding="utf-8").read().splitlines()[0][2:].strip()
    p = doc.add_paragraph()
    r = p.add_run(title); r.bold = True; r.font.size = Pt(15); r.font.name = FONT
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    a = doc.add_paragraph()
    ra = a.add_run("Yongxin Yang, B.M."); ra.font.size = Pt(12); ra.font.name = FONT
    a.alignment = WD_ALIGN_PARAGRAPH.CENTER
    aff = doc.add_paragraph()
    rf = aff.add_run("The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China")
    rf.font.size = Pt(11); rf.font.name = FONT; rf.italic = True
    aff.alignment = WD_ALIGN_PARAGRAPH.CENTER
    corr = doc.add_paragraph()
    rc = corr.add_run("*Corresponding author. Email: 960856791@qq.com; ORCID: 0009-0004-9698-6552*")
    rc.font.size = Pt(11); rc.font.name = FONT
    corr.alignment = WD_ALIGN_PARAGRAPH.CENTER
    sh = doc.add_paragraph()
    rs = sh.add_run("Short title: Opioid transcriptional silence in spinal microglia")
    rs.font.size = Pt(10); rs.font.name = FONT; rs.italic = True
    sh.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_paragraph()

# ---------- manuscript conversion ----------
def in_author_preamble(line_idx, lines):
    # author-note preamble sits between the H1 title and the first '---' + Abstract
    return False  # handled by callers

def build_manuscript():
    raw = io.open(SRC, "r", encoding="utf-8").read()
    # strip HTML-style nothing; split lines
    lines = raw.splitlines()
    doc = Document()
    set_base_style(doc)
    inject_title_page(doc)

    i = 0
    n = len(lines)
    skip_author_note = False  # True while inside §6 "Author's note on journals" block
    started = False  # only convert from the Abstract section onward
    table_buf = []

    def flush_table(buf):
        if len(buf) < 2:
            return
        # buf[0] = header, buf[1] = separator, rest = rows
        header = [c.strip() for c in buf[0].strip("|").split("|")]
        rows = []
        for bl in buf[2:]:
            rows.append([c.strip() for c in bl.strip("|").split("|")])
        t = doc.add_table(rows=1, cols=len(header))
        t.style = "Light Grid Accent 1"
        hdr = t.rows[0].cells
        for j, htext in enumerate(header):
            hdr[j].text = htext
            for p in hdr[j].paragraphs:
                for r in p.runs:
                    r.font.name = FONT; r.font.size = Pt(10); r.bold = True
        for row in rows:
            cells = t.add_row().cells
            for j, ctext in enumerate(row):
                cells[j].text = ctext
                for p in cells[j].paragraphs:
                    for r in p.runs:
                        r.font.name = FONT; r.font.size = Pt(10)
        doc.add_paragraph()

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # detect start: the Abstract section
        if not started:
            if stripped.startswith("## ") and "Abstract" in stripped:
                started = True
            else:
                i += 1
                continue

        # toggle §6 author-note skip
        if stripped.startswith("> **Author's note on journals"):
            skip_author_note = True
            i += 1
            continue
        if skip_author_note:
            if stripped.startswith(">"):
                i += 1
                continue
            else:
                skip_author_note = False  # block ended

        if skip_author_note:
            i += 1
            continue

        # blank
        if stripped == "":
            if table_buf:
                flush_table(table_buf); table_buf = []
            i += 1
            continue

        # horizontal rule
        if stripped == "---":
            if table_buf:
                flush_table(table_buf); table_buf = []
            i += 1
            continue

        # table row
        if stripped.startswith("|") and stripped.endswith("|"):
            # could be separator
            if set(stripped.replace("|", "").replace("-", "").replace(":", "").replace(" ", "")) == set():
                # separator line; keep header in buffer if present
                if table_buf:
                    table_buf.append(stripped)  # will be ignored by flush (needs>=2 and uses buf[2:])
                i += 1
                continue
            if table_buf:
                flush_table(table_buf); table_buf = []
            table_buf.append(stripped)
            i += 1
            continue
        else:
            if table_buf:
                flush_table(table_buf); table_buf = []

        # headings
        if stripped.startswith("# "):
            add_heading(doc, stripped[2:].strip(), 1); i += 1; continue
        if stripped.startswith("## "):
            add_heading(doc, stripped[3:].strip(), 1); i += 1; continue
        if stripped.startswith("### "):
            add_heading(doc, stripped[4:].strip(), 2); i += 1; continue

        # blockquote (manuscript callout, e.g. "Reversal of the original hypothesis")
        if stripped.startswith(">"):
            txt = stripped.lstrip("> ").strip()
            if txt == "":
                i += 1; continue
            add_para(doc, txt, italic=True, size=11, color=RGBColor(0x55, 0x55, 0x55))
            i += 1; continue

        # bullet / numbered list
        m = re.match(r"^(\d+)\.\s+(.*)$", stripped)
        if m:
            p = doc.add_paragraph(style="List Number")
            r = p.add_run(m.group(2)); r.font.name = FONT; r.font.size = Pt(12)
            i += 1; continue
        if stripped.startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            r = p.add_run(stripped[2:].strip()); r.font.name = FONT; r.font.size = Pt(12)
            i += 1; continue

        # plain paragraph
        add_para(doc, stripped, size=12)
        i += 1

    if table_buf:
        flush_table(table_buf)

    # figures: copy PNGs to figures/ and add a reference list
    pngs = sorted(glob.glob(os.path.join(FIGS, "Fig*.png")))
    if pngs:
        doc.add_page_break()
        add_heading(doc, "Figures", 1)
        for p in pngs:
            base = os.path.basename(p)
            shutil.copy(p, os.path.join(FIGOUT, base))
            cap = "Figure %s — see manuscript §3 for description." % base.replace("Fig", "").replace("_decoy_null", " (decoy null)").replace(".png", "")
            add_para(doc, cap, italic=True, size=10, color=RGBColor(0x55,0x55,0x55))

    # strip residual markdown markers from runs
    for p in doc.paragraphs:
        for r in p.runs:
            r.text = r.text.replace("**", "").replace("`", "").replace("__", "")
    for t in doc.tables:
        for row in t.rows:
            for c in row.cells:
                for p in c.paragraphs:
                    for r in p.runs:
                        r.text = r.text.replace("**", "").replace("`", "").replace("__", "")

    out_docx = os.path.join(OUT, "manuscript_plosone_v5.docx")
    doc.save(out_docx)
    print("[build_submission] wrote", out_docx)

# ---------- cover letter (PLOS ONE: <=1 page, no APC mention) ----------
def build_cover_letter():
    md = """# Cover Letter

**To:** PLOS ONE Editorial Office (journals.plos.org/plosone)

**Article type:** Research Article

**Title:** Opioid exposure is transcriptionally silent in spinal cord microglia: a multi-omics dissection of the neuropathic-pain activation hub universe and its RUNX1-anchored, repurposable-analgesia target space

---

Dear PLOS ONE Editors,

We submit the above Research Article for consideration in *PLOS ONE*.

**Contribution to the literature.** We provide a rigorously gated, reproducible re-analysis of public mouse spinal-cord microglia transcriptomics (GSE117320) showing that long-term opioid (morphine) exposure is transcriptionally silent in microglia, whereas neuropathic pain drives a 3,069-gene activation program converged on six hub genes and mechanistically anchored to RUNX1 and the NF-kappaB axis. We add honestly executed validation layers (a LINCS positive-control recovery gate and structure-based AutoDock Vina docking with a decoy/random-box specificity panel) and nominate phenazone as a brain-penetrant repurposing lead. The work is notable for explicitly reporting negative and modest results (G1 NOT MET; docking specificity p~0.02) rather than overstating them.

**Relation to prior work.** This is an original re-analysis of publicly deposited data; it is not under consideration elsewhere and reports no prior duplicate publication. The conceptual frame extends our group's multi-omics virtual-kockout pipeline to opioid tolerance / opioid-induced hyperalgesia.

**Why PLOS ONE.** *PLOS ONE* evaluates submissions on scientific soundness, not perceived novelty or impact. Our manuscript is a methodologically rigorous computational study that openly reports honest negatives and modest positives; it is squarely within *PLOS ONE*'s scope (neuroscience, genomics, computational biology & bioinformatics, pharmacology).

**Prior interaction with PLOS.** None.

**Suggested Academic Editors.** [Suggest 2-3 editors from the PLOS ONE board with expertise in neuroinflammation / glial biology / computational drug repurposing; fill before submission.]

**Opposed reviewers.** None.

We confirm that all underlying data are publicly available (GEO accessions GSE117320 and GSE92742) and that a Data Availability Statement, Author Contributions (CRediT), Competing Interests, and Funding statements are included in the manuscript.

Thank you for your consideration.

Sincerely,

Yongxin Yang, B.M.
The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine
Fuzhou, Fujian 350003, China
960856791@qq.com
"""
    out_md = os.path.join(OUT, "cover_letter.md")
    io.open(out_md, "w", encoding="utf-8").write(md)
    # build docx
    doc = Document(); set_base_style(doc)
    for ln in md.splitlines():
        s = ln.strip()
        if s == "": continue
        if s.startswith("# "):
            add_heading(doc, s[2:].strip(), 1); continue
        if s == "---": continue
        if s.startswith("**") and s.endswith("**"):
            add_para(doc, s.strip("*"), bold=True, size=12); continue
        add_para(doc, s, size=12)
    out_docx = os.path.join(OUT, "cover_letter.docx")
    doc.save(out_docx)
    print("[build_submission] wrote", out_docx, "and", out_md)

# ---------- highlights ----------
def build_highlights():
    items = [
        "Opioid (morphine) exposure is transcriptionally silent in mouse spinal microglia.",
        "Neuropathic pain drives a 3,069-gene microglial activation hub universe (6 hubs).",
        "The hub universe is anchored to RUNX1 (ChEA_2016 p=3.4x10-4, OR=1.57).",
        "Reversal landscape converges on NF-kB axis; phenazone is a brain-penetrant lead.",
        "Honest reporting: G1/G4 NOT MET; docking specific but modest (decoy p~0.02).",
    ]
    md = "# Highlights\n\n" + "\n".join("- " + it for it in items) + "\n"
    out_md = os.path.join(OUT, "highlights.md")
    io.open(out_md, "w", encoding="utf-8").write(md)
    doc = Document(); set_base_style(doc)
    add_heading(doc, "Highlights", 1)
    for it in items:
        p = doc.add_paragraph(style="List Bullet")
        r = p.add_run(it); r.font.name = FONT; r.font.size = Pt(12)
    out_docx = os.path.join(OUT, "highlights.docx")
    doc.save(out_docx)
    print("[build_submission] wrote", out_docx, "and", out_md)

# ---------- response letter template ----------
def build_response_letter():
    md = """# Response to Reviewers (template)

**Manuscript:** Opioid exposure is transcriptionally silent in spinal cord microglia...
**Reference ID:** [journal ms ID]

Dear Editors and Reviewers,

We thank the reviewers for their constructive comments. We have addressed each point below (our responses in plain text; manuscript changes noted by section).

---

## Reviewer 1

**Comment 1.** [paste reviewer comment]

**Response.** [your response — be specific; cite the revised line/section and any new number]

**Change.** [e.g., §3.10 revised to report decoy/random-box p-values; Fig9 added]

---

## Reviewer 2

**Comment 1.** [paste reviewer comment]

**Response.** [your response]

**Change.** [note]

---

## Summary of manuscript changes
- [list files / sections revised]
- [note any new analyses added (e.g., decoy specificity panel)]
"""
    out_md = os.path.join(OUT, "response_letter_template.md")
    io.open(out_md, "w", encoding="utf-8").write(md)
    doc = Document(); set_base_style(doc)
    for ln in md.splitlines():
        s = ln.strip()
        if s == "": continue
        if s.startswith("# "):
            add_heading(doc, s[2:].strip(), 1); continue
        if s == "---": continue
        if s.startswith("**") and s.endswith("**"):
            add_para(doc, s.strip("*"), bold=True, size=12); continue
        add_para(doc, s, size=12)
    out_docx = os.path.join(OUT, "response_letter_template.docx")
    doc.save(out_docx)
    print("[build_submission] wrote", out_docx, "and", out_md)

# ---------- reporting checklist / transparency statement ----------
def build_reporting_checklist():
    md = """# Reporting checklist & transparency statement (PLOS ONE submission)

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
"""
    out_md = os.path.join(OUT, "reporting_checklist.md")
    io.open(out_md, "w", encoding="utf-8").write(md)
    doc = Document(); set_base_style(doc)
    for ln in md.splitlines():
        s = ln.strip()
        if s == "": continue
        if s.startswith("# "):
            add_heading(doc, s[2:].strip(), 1); continue
        if s.startswith("**") and s.endswith("**"):
            add_para(doc, s.strip("*"), bold=True, size=12); continue
        add_para(doc, s, size=12)
    out_docx = os.path.join(OUT, "reporting_checklist.docx")
    doc.save(out_docx)
    print("[build_submission] wrote", out_docx, "and", out_md)

if __name__ == "__main__":
    build_manuscript()
    build_cover_letter()
    build_highlights()
    build_response_letter()
    build_reporting_checklist()
    print("[build_submission] DONE")
