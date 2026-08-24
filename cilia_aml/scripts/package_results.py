#!/usr/bin/env python3
"""Assemble every deliverable into one organised, zipped bundle in Downloads."""
import os
import shutil
import zipfile

DL = "/mnt/c/Users/imabb/Downloads"
OUT = os.path.join(DL, "Transcriptomics_Results_Complete")
shutil.rmtree(OUT, ignore_errors=True)

# (source, destination-relative-path)
ITEMS = [
    # ---------- AML / miR-551b : reports ----------
    ("cilia_aml/reports/Cilia_AML_miR551b_report.pdf",  "01_AML_miR551b/reports/Cilia_AML_miR551b_report.pdf"),
    ("cilia_aml/reports/Cilia_AML_miR551b_report.docx", "01_AML_miR551b/reports/Cilia_AML_miR551b_report.docx"),
    ("cilia_aml/reports/Cilia_AML_miR551b_report.html", "01_AML_miR551b/reports/Cilia_AML_miR551b_report.html"),
    # ---------- AML : HD figures ----------
    ("cilia_aml/figures/hd/pca.png",                 "01_AML_miR551b/figures_300dpi/AML_01_PCA.png"),
    ("cilia_aml/figures/hd/volcano.png",             "01_AML_miR551b/figures_300dpi/AML_02_volcano_genomewide.png"),
    ("cilia_aml/figures/hd/ma.png",                  "01_AML_miR551b/figures_300dpi/AML_03_MA_plot.png"),
    ("cilia_aml/figures/hd/volcano_cilia_hd.png",    "01_AML_miR551b/figures_300dpi/AML_04_volcano_cilia_683set.png"),
    ("cilia_aml/figures/hd/volcano_network_hd.png",  "01_AML_miR551b/figures_300dpi/AML_05_volcano_network_25genes.png"),
    # ---------- AML : tables ----------
    ("cilia_aml/de_paired/de_results.csv",        "01_AML_miR551b/tables/AML_DE_paired_model_ALL_GENES.csv"),
    ("cilia_aml/de_primary/de_results.csv",       "01_AML_miR551b/tables/AML_DE_primaryAML_only.csv"),
    ("cilia_aml/de_kg1a/de_results.csv",          "01_AML_miR551b/tables/AML_DE_KG1a_only.csv"),
    ("cilia_aml/de_paired/cilia_in_aml.xlsx",     "01_AML_miR551b/tables/cilia_683set_in_AML.xlsx"),
    ("cilia_aml/de_paired/network_in_aml.xlsx",   "01_AML_miR551b/tables/network_25set_in_AML.xlsx"),
    ("cilia_aml/de_paired/IFT_BBSome_in_AML.xlsx", "01_AML_miR551b/tables/IFT_BBSome_in_AML.xlsx"),
    ("cilia_aml/de_paired/IFT_BBSome_in_AML.csv",  "01_AML_miR551b/tables/IFT_BBSome_in_AML.csv"),
    ("cilia_aml/data/counts_gene.tsv",            "01_AML_miR551b/tables/AML_gene_count_matrix.tsv"),
    ("cilia_aml/data/design.csv",                 "01_AML_miR551b/tables/AML_design.csv"),
    # ---------- human GSE336399 ----------
    ("results/hd/gse336399/pca.png",                    "02_GSE336399_human_nasal/figures_300dpi/HUMAN_01_PCA.png"),
    ("results/hd/gse336399/volcano.png",                "02_GSE336399_human_nasal/figures_300dpi/HUMAN_02_volcano_genomewide.png"),
    ("results/hd/gse336399/ma.png",                     "02_GSE336399_human_nasal/figures_300dpi/HUMAN_03_MA_plot.png"),
    ("results/hd/gse336399/volcano_cilia_hd.png",       "02_GSE336399_human_nasal/figures_300dpi/HUMAN_04_volcano_cilia_683set.png"),
    ("results/hd/gse336399/volcano_network_hd.png",     "02_GSE336399_human_nasal/figures_300dpi/HUMAN_05_volcano_network_25genes.png"),
    ("results/figures/GSE336399_human_enrichment.png",  "02_GSE336399_human_nasal/figures_300dpi/HUMAN_06_enrichment.png"),
    ("results/gse336399/de_results.csv",                "02_GSE336399_human_nasal/tables/HUMAN_DE_results.csv"),
    ("results/gse336399/GSE336399_results.xlsx",        "02_GSE336399_human_nasal/tables/GSE336399_results.xlsx"),
    ("results/gse336399/GSE336399_geneset_hits.xlsx",   "02_GSE336399_human_nasal/tables/GSE336399_cilia_geneset_hits.xlsx"),
    ("results/gse336399/enrichr_up.csv",                "02_GSE336399_human_nasal/tables/enrichment_up.csv"),
    ("results/gse336399/enrichr_down.csv",              "02_GSE336399_human_nasal/tables/enrichment_down.csv"),
    # ---------- mouse GSE334363 ----------
    ("results/hd/gse334363/pca.png",                    "03_GSE334363_mouse_sildenafil/figures_300dpi/MOUSE_01_PCA.png"),
    ("results/hd/gse334363/volcano.png",                "03_GSE334363_mouse_sildenafil/figures_300dpi/MOUSE_02_volcano_genomewide.png"),
    ("results/hd/gse334363/ma.png",                     "03_GSE334363_mouse_sildenafil/figures_300dpi/MOUSE_03_MA_plot.png"),
    ("results/hd/gse334363/volcano_cilia_hd.png",       "03_GSE334363_mouse_sildenafil/figures_300dpi/MOUSE_04_volcano_cilia_683set.png"),
    ("results/hd/gse334363/volcano_network_hd.png",     "03_GSE334363_mouse_sildenafil/figures_300dpi/MOUSE_05_volcano_network_25genes.png"),
    ("results/figures/GSE334363_mouse_enrichment.png",  "03_GSE334363_mouse_sildenafil/figures_300dpi/MOUSE_06_enrichment.png"),
    ("results/gse334363/de_results.csv",                "03_GSE334363_mouse_sildenafil/tables/MOUSE_DE_results.csv"),
    ("results/gse334363/GSE334363_results.xlsx",        "03_GSE334363_mouse_sildenafil/tables/GSE334363_results.xlsx"),
    ("results/gse334363/GSE334363_geneset_hits.xlsx",   "03_GSE334363_mouse_sildenafil/tables/GSE334363_cilia_geneset_hits.xlsx"),
    ("results/gse334363/enrichr_up.csv",                "03_GSE334363_mouse_sildenafil/tables/enrichment_up.csv"),
    ("results/gse334363/enrichr_down.csv",              "03_GSE334363_mouse_sildenafil/tables/enrichment_down.csv"),
    # ---------- cross-dataset + earlier report ----------
    ("cilia_aml/de/network_genes_in_existing.csv", "04_cross_dataset/network_25genes_human_vs_mouse.csv"),
    ("results/geneset_in_both.xlsx",               "04_cross_dataset/cilia_683set_human_vs_mouse.xlsx"),
    ("results/Ciliary_geneset_report.pdf",         "04_cross_dataset/Ciliary_geneset_report_2datasets.pdf"),
    ("results/Ciliary_geneset_report.docx",        "04_cross_dataset/Ciliary_geneset_report_2datasets.docx"),
    ("results/Ciliary_geneset_report.html",        "04_cross_dataset/Ciliary_geneset_report_2datasets.html"),
    # ---------- gene sets used ----------
    ("cilia_aml/data/cilia_genes.txt",   "05_gene_sets/cilia_683_genes.txt"),
    ("cilia_aml/data/network_genes.txt", "05_gene_sets/interaction_network_25_genes.txt"),
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

README = """TRANSCRIPTOMICS RESULTS - COMPLETE BUNDLE
=========================================
Generated by the open transcriptomics pipeline (github.com/imabbi47/transcriptomics)

01_AML_miR551b/          NEW: anti-miR-551b vs Control in AML (BioProject PRJNA1513000)
                         Processed end-to-end from RAW READS: aria2 download from the SRA
                         AWS mirror -> Salmon selective alignment vs GENCODE v44
                         (72-83% mapping) -> DESeq2. 8 single-end runs, 2x2x2 design
                         (primary AML / KG1a  x  LNA / ZIP  x  anti-miR / control).
    RESULT: near-null. 3 significant genes genome-wide (TCF7, IGKC, 1 unannotated).
            Cilia 683-set: 0/563 significant. Network 25-set: 0/17 significant.
            IFT/BBSome: 41/46 expressed, 0 changed (not even nominal p<0.05).
            Per-cell-type: primary AML gave 33 hits but dominated by T/B-cell
            markers = immune-composition confound, not miR-551b. KG1a ~0.

02_GSE336399_human_nasal/   Human nasal brushings, Cancer vs Benign (92 samples).
    RESULT: 479 DE genes (16 up / 463 down). Inflammatory + interferon suppression.
            Cilia set coordinately DOWN: 8 significant, all down-regulated.

03_GSE334363_mouse_sildenafil/  Mouse tumour lines, Sildenafil vs Control (18 samples).
    RESULT: 6,553 DE genes (3,397 up / 3,156 down). Cholesterol homeostasis +
            interferon UP; Myc targets + ribosome biogenesis DOWN.
            Cilia set: 223 significant but bidirectional (128 up / 95 down) = incidental.

04_cross_dataset/        Cilia set and interaction network compared across datasets,
                         plus the earlier two-dataset journal report (PDF/DOCX/HTML).

05_gene_sets/            The exact gene lists used (683 cilia genes; 25 network genes).

NOTES
-----
* All figures are 300 dpi (1980x1500 to 2280x2400 px), publication-ready raster.
* Reports contain full Introduction / Methods / Results / Discussion + limitations.
* The human and mouse analyses used the authors' processed counts (raw reads were
  embargoed); the AML analysis used raw reads. Human/mouse results reproduced
  byte-identically on re-run.
* KEY CAVEAT for AML: the design gives one biological unit per unique condition
  (LNA/ZIP chemistries act as pseudo-replicates), so power to detect modest
  effects is low. The negative result should be read with that in mind.
"""
open(os.path.join(OUT, "README.txt"), "w", encoding="utf-8").write(README)

# zip it
zip_path = os.path.join(DL, "Transcriptomics_Results_Complete.zip")
if os.path.exists(zip_path):
    os.remove(zip_path)
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    for root, _dirs, files in os.walk(OUT):
        for f in files:
            p = os.path.join(root, f)
            z.write(p, os.path.relpath(p, OUT))

print("copied {} files".format(copied))
if missing:
    print("MISSING ({}):".format(len(missing)))
    for m in missing:
        print("   " + m)
print("\nbundle : {}".format(OUT))
print("zip    : {}  ({:.1f} MB)".format(zip_path, os.path.getsize(zip_path) / 1e6))
print("\nstructure:")
for root, _dirs, files in os.walk(OUT):
    rel = os.path.relpath(root, OUT)
    depth = 0 if rel == "." else rel.count(os.sep) + 1
    name = "." if rel == "." else os.path.basename(root)
    print("  " * depth + name + "/")
    for f in sorted(files):
        print("  " * (depth + 1) + f)
