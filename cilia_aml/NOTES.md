# Cilia in AML - transcriptomics

Goal: extend the cilia / ciliopathy gene-set analysis to Acute Myeloid
Leukemia (AML) - testing whether ciliary and centrosomal genes are
dysregulated in AML versus normal hematopoiesis.

Note: mature blood cells are largely non-ciliated, so in AML the signal is
more likely in centrosome / basal-body genes than motile-cilia genes -
a specific, less-obvious angle worth checking carefully.

Candidate data sources:
- TCGA-LAML (via UCSC Xena or recount3) - AML tumors
- GEO AML RNA-seq series (AML blasts vs normal CD34+ / bone marrow)
- BeatAML cohort

Folder layout:
  data/     count matrices, designs, references
  de/       differential-expression outputs
  figures/  PCA, volcano, enrichment, heatmaps
  reports/  HTML / DOCX / PDF write-ups
