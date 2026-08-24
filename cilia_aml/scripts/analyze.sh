#!/usr/bin/env bash
# Downstream analysis for Cilia+AML once counts_gene.tsv exists:
#   DESeq2 -> cilia & network cross-refs -> highlighted volcanos -> enrichment
set +e
PROJ=/home/abhi/github-projects/brainstroming_transcriptomics
cd "$PROJ" || exit 1
PY="$PROJ/.venv/bin/python"
D=cilia_aml/data
DE=cilia_aml/de
FIG=cilia_aml/figures
LOG=$D/analyze.log
: > "$LOG"
say(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
mkdir -p "$DE" "$FIG"

if [ ! -s "$D/counts_gene.tsv" ]; then say "FATAL: counts_gene.tsv missing"; exit 1; fi

say "=== DESeq2: miR-551b-inhibited vs Control (adjusted for cell type), n=4 vs 4 ==="
$PY transcriptomics/de.py --counts "$D/counts_gene.tsv" --design "$D/design.csv" \
  --factor treatment --ref Control --alt Mir551b --covariate cell_type \
  --outdir "$DE" >> "$LOG" 2>&1
say "de rc=$?"
cp -f "$DE/volcano.png" "$FIG/volcano_all.png" 2>/dev/null
cp -f "$DE/pca.png" "$FIG/pca.png" 2>/dev/null
cp -f "$DE/ma.png" "$FIG/ma.png" 2>/dev/null

say "=== cilia 683-gene set cross-ref ==="
$PY analysis/lookup_geneset.py --genes "$D/cilia_genes.txt" \
  --de "AML_miR551bKD:$DE/de_results.csv" --out "$DE/cilia_in_aml.xlsx" >> "$LOG" 2>&1
say "cilia lookup rc=$?"

say "=== network 25-gene set cross-ref ==="
$PY analysis/lookup_geneset.py --genes "$D/network_genes.txt" \
  --de "AML_miR551bKD:$DE/de_results.csv" --out "$DE/network_in_aml.xlsx" >> "$LOG" 2>&1
say "network lookup rc=$?"

say "=== volcano: cilia set highlighted ==="
$PY analysis/volcano_geneset.py --de "$DE/de_results.csv" --genes "$D/cilia_genes.txt" \
  --title "Cilia/centrosome genes: miR-551b inhibited vs Control (AML+KG1a)" \
  --out "$FIG/volcano_cilia.png" >> "$LOG" 2>&1
say "=== volcano: network set highlighted ==="
$PY analysis/volcano_geneset.py --de "$DE/de_results.csv" --genes "$D/network_genes.txt" \
  --title "Interaction-network genes: miR-551b inhibited vs Control (AML+KG1a)" \
  --out "$FIG/volcano_network.png" >> "$LOG" 2>&1

say "=== enrichment (Enrichr: Hallmark / KEGG / GO-BP) ==="
$PY transcriptomics/enrich.py --results "$DE/de_results.csv" --outdir "$DE" \
  --label "miR-551b inhibited vs Control" >> "$LOG" 2>&1
say "enrich rc=$?"
cp -f "$DE"/enrichment_up.png "$FIG/" 2>/dev/null
cp -f "$DE"/enrichment_down.png "$FIG/" 2>/dev/null

say "ANALYZE DONE"
