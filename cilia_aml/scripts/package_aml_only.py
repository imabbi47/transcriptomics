#!/usr/bin/env python3
"""Assemble the AML (miR-551b) results only, into an organised zipped bundle."""
import json
import os
import shutil
import zipfile

import pandas as pd

DL = "/mnt/c/Users/imabb/Downloads"
OUT = os.path.join(DL, "AML_miR551b_Transcriptomics")
shutil.rmtree(OUT, ignore_errors=True)

RUNS = {
    "SRR40166669": ("primary_AML", "LNA", "Mir551b"),
    "SRR40166670": ("primary_AML", "LNA", "Control"),
    "SRR40166671": ("primary_AML", "ZIP", "Mir551b"),
    "SRR40166672": ("primary_AML", "ZIP", "Control"),
    "SRR40166673": ("KG1a", "LNA", "Mir551b"),
    "SRR40166674": ("KG1a", "LNA", "Control"),
    "SRR40166675": ("KG1a", "ZIP", "Mir551b"),
    "SRR40166676": ("KG1a", "ZIP", "Control"),
}

# ---------- QC table straight from the Salmon run metadata ----------
qc = []
for run, (cell, chem, trt) in RUNS.items():
    meta = "cilia_aml/data/quant/{}/aux_info/meta_info.json".format(run)
    row = dict(run=run, cell_type=cell, chemistry=chem, treatment=trt)
    if os.path.exists(meta):
        m = json.load(open(meta))
        row["reads_processed"] = m.get("num_processed")
        row["reads_mapped"] = m.get("num_mapped")
        rate = m.get("percent_mapped")
        row["mapping_rate_pct"] = round(rate, 2) if isinstance(rate, (int, float)) else None
        row["library_type"] = m.get("library_types", [None])[0]
    qc.append(row)
qc = pd.DataFrame(qc)
os.makedirs(os.path.join(OUT, "tables"), exist_ok=True)
qc.to_csv(os.path.join(OUT, "tables", "AML_sequencing_QC.csv"), index=False)

ITEMS = [
    # reports
    ("cilia_aml/reports/Cilia_AML_miR551b_report.pdf",  "reports/AML_miR551b_report.pdf"),
    ("cilia_aml/reports/Cilia_AML_miR551b_report.docx", "reports/AML_miR551b_report.docx"),
    ("cilia_aml/reports/Cilia_AML_miR551b_report.html", "reports/AML_miR551b_report.html"),
    # figures (300 dpi)
    ("cilia_aml/figures/hd/pca.png",                "figures_300dpi/AML_01_PCA.png"),
    ("cilia_aml/figures/hd/volcano.png",            "figures_300dpi/AML_02_volcano_genomewide.png"),
    ("cilia_aml/figures/hd/ma.png",                 "figures_300dpi/AML_03_MA_plot.png"),
    ("cilia_aml/figures/hd/volcano_cilia_hd.png",   "figures_300dpi/AML_04_volcano_cilia_683set.png"),
    ("cilia_aml/figures/hd/volcano_network_hd.png", "figures_300dpi/AML_05_volcano_network_25genes.png"),
    # DE tables
    ("cilia_aml/de_paired/de_results.csv",  "tables/AML_DE_paired_model_ALL_GENES.csv"),
    ("cilia_aml/de_primary/de_results.csv", "tables/AML_DE_primaryAML_only.csv"),
    ("cilia_aml/de_kg1a/de_results.csv",    "tables/AML_DE_KG1a_only.csv"),
    # gene-set results
    ("cilia_aml/de_paired/cilia_in_aml.xlsx",     "tables/cilia_683set_in_AML.xlsx"),
    ("cilia_aml/de_paired/network_in_aml.xlsx",   "tables/network_25set_in_AML.xlsx"),
    ("cilia_aml/de_paired/IFT_BBSome_in_AML.xlsx", "tables/IFT_BBSome_in_AML.xlsx"),
    ("cilia_aml/de_paired/IFT_BBSome_in_AML.csv",  "tables/IFT_BBSome_in_AML.csv"),
    # inputs / provenance
    ("cilia_aml/data/counts_gene.tsv", "tables/AML_gene_count_matrix.tsv"),
    ("cilia_aml/data/design.csv",      "tables/AML_design.csv"),
    ("cilia_aml/data/cilia_genes.txt",   "gene_sets/cilia_683_genes.txt"),
    ("cilia_aml/data/network_genes.txt", "gene_sets/interaction_network_25_genes.txt"),
    ("cilia_aml/data/pipeline.log", "provenance/pipeline_run.log"),
    ("cilia_aml/data/analyze.log",  "provenance/analysis_run.log"),
    ("cilia_aml/scripts/pipeline.sh",       "provenance/scripts/pipeline.sh"),
    ("cilia_aml/scripts/aggregate.py",      "provenance/scripts/aggregate.py"),
    ("cilia_aml/scripts/analyze.sh",        "provenance/scripts/analyze.sh"),
    ("cilia_aml/scripts/ift_bbsome.py",     "provenance/scripts/ift_bbsome.py"),
    ("cilia_aml/scripts/make_aml_report.py", "provenance/scripts/make_aml_report.py"),
]

copied, missing = 0, []
for src, rel in ITEMS:
    dst = os.path.join(OUT, rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(src):
        shutil.copy2(src, dst)
        copied += 1
    else:
        missing.append(src)

README = """AML / miR-551b TRANSCRIPTOMICS - RESULTS
=======================================
Dataset : NCBI BioProject PRJNA1513000, runs SRR40166669-SRR40166676
Design  : 2x2x2 - cell context (primary AML blasts | KG1a cell line)
                x antisense chemistry (LNA | ZIP)
                x treatment (anti-miR-551b | control oligo)
Question: does inhibiting miR-551b change the AML transcriptome, and in
          particular the ciliary/centrosomal gene programme?

PROCESSING (from raw reads)
---------------------------
Download    : aria2c, 16 connections, SRA AWS Open-Data mirror (full depth,
              13-25 M single-end reads/run; ~4 GB total)
Quantify    : Salmon 2.5.1 selective alignment vs GENCODE v44 (251,955 transcripts)
Mapping rate: 71.8-72.7% (primary AML), 80.7-82.9% (KG1a) - see tables/AML_sequencing_QC.csv
Matrix      : 60,883 genes x 8 samples; 25,994 passed total-count>=10 and were tested
Statistics  : DESeq2 (pyDESeq2), Wald test, Benjamini-Hochberg FDR
Primary model: ~ cell_type + chemistry + treatment  (paired/blocked)

RESULTS
-------
Genome-wide (paired model): 3 significant genes at FDR<0.05, all DOWN
    ENSG00000257379  -8.46  (FDR 1.6e-3)   unannotated locus
    TCF7             -2.16  (FDR 1.7e-2)   Wnt / T-cell transcription factor
    IGKC             -1.71  (FDR 1.5e-2)   immunoglobulin kappa constant

Cilia 683-gene set      : 563 expressed, 0 significant.
                          Mean |log2FC| 0.185 vs 0.321 background -> MORE stable
                          than the average gene.
                          Motile/axonemal genes sit at the noise floor (median 14
                          counts); the expressed cilia programme in AML is
                          centrosomal/basal-body (median 451) and trafficking (527).

IFT / BBSome machinery  : 41/46 expressed, 0 significant - not even one gene at
                          nominal p<0.05 (you would expect ~2 by chance).
                          Mean |log2FC| 0.104, i.e. ~3x more stable than average.
                          Largest shift BBS1 +0.52 (p=0.093, FDR 1.00).
                          IFT88 -0.04, ARL13B +0.21 - both flat.
                          NOT detected: WDR34, WDR60, TCTEX1D2 (dynein-2 accessory
                          subunits), HSPB11, RABL4 -> the retrograde motor is
                          incomplete, as expected for a non-ciliated cell.

Interaction network (25): 17 expressed, 0 significant. The 8 non-detected members
                          were the dopaminergic module (SLC6A3, DBH, ANKK1, DRD1/2/5,
                          PPP1R1B) plus IL13 - neuronal genes, not expressed in myeloid
                          cells. TLR4 +0.47 was the largest move (ns).

Per cell type (secondary):
    primary AML (2v2): 33 significant genes - BUT led by T/B-lymphocyte markers
                       (TRBC1, ZAP70, LCK, IGKC, IGHM, RORA, TCF7). With n=1 per
                       unique condition and variable normal immune-cell content in
                       patient material, this is best explained as a cell-composition
                       difference between samples, NOT a miR-551b programme.
    KG1a (2v2)       : 3 hits, dominated by read-through/unannotated loci -> effectively
                       no reproducible effect in the clean cell line.

CONCLUSION
----------
miR-551b inhibition produces a minimal, non-reproducible transcriptional footprint in
AML cells and does not engage the ciliary/centrosomal programme. The cilia machinery is
present and transcribed in AML, but it is transcriptionally inert here and structurally
incomplete (missing dynein-2 accessory subunits). This is a NEGATIVE result.

KEY LIMITATION
--------------
The design provides only ONE biological unit per unique condition; the LNA and ZIP
chemistries function as pseudo-replicates. Power to detect modest coordinated effects is
therefore low, and a subtle shift could hide below detection. Additionally, miRNA
inhibition often acts translationally, so a small mRNA footprint does not exclude a
protein-level effect. Primary-AML samples carry variable normal-cell content.

FOLDER GUIDE
------------
reports/          Full report (Intro/Methods/Results/Discussion + limitations) - PDF, DOCX, HTML
figures_300dpi/   PCA, genome-wide volcano, MA plot, cilia-set volcano, network volcano
tables/           DE results (all genes + per cell type), gene-set results, IFT/BBSome,
                  the gene count matrix, design, and sequencing QC
gene_sets/        The exact gene lists used
provenance/       Run logs and the actual scripts used, for reproducibility
"""
open(os.path.join(OUT, "README.txt"), "w", encoding="utf-8").write(README)

zip_path = os.path.join(DL, "AML_miR551b_Transcriptomics.zip")
if os.path.exists(zip_path):
    os.remove(zip_path)
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    for root, _d, files in os.walk(OUT):
        for f in files:
            p = os.path.join(root, f)
            z.write(p, os.path.relpath(p, OUT))

print("copied {} files (+ QC table + README)".format(copied))
if missing:
    print("MISSING:")
    for m in missing:
        print("   " + m)
print("\nSEQUENCING QC")
print(qc.to_string(index=False))
print("\nbundle : {}".format(OUT))
print("zip    : {} ({:.1f} MB)".format(zip_path, os.path.getsize(zip_path) / 1e6))
