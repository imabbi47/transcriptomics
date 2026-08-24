# Project status — transcriptomics pipeline

**Repo:** github.com/imabbi47/transcriptomics (5 commits, `main` in sync)
**Local:** `/home/abhi/github-projects/brainstroming_transcriptomics` (WSL Ubuntu)

## Built (all committed) — CLI `transcriptomics <cmd>`
- **resolve / fetch** — accession (GEO/SRA/ENA/DDBJ) → runs → FASTQ + samplesheet (stdlib core)
- **quantify** — Salmon quasi-mapping (alignment-free)
- **align** — HISAT2/STAR → featureCounts (alignment-based)
- **de** — differential expression (pyDESeq2)
- **enrich** — GO/KEGG/Hallmark (gseapy/Enrichr)
- **report / make_docs** — single-file HTML, plus HTML/DOCX/PDF
- **geo-design** — parse GEO series_matrix → design table
- **serve** — FastAPI web backend (`webapp.py`) + static site (`docs/`, live at transcriptomics.vercel.app)
- helpers in `analysis/`: export_excel, lookup_geneset, volcano_geneset, make_figures, prep_gse336399

## Datasets analysed (results in `results/` [gitignored] + `Downloads/Transcriptomics_results/`)
- **GSE334363** (mouse, Sildenafil vs Control): 6,553 DE; cholesterol/interferon ↑, Myc/ribosome ↓
- **GSE336399** (human nasal, Cancer vs Benign, 92 samples): 479 DE (463 ↓); immune/interferon suppression
- **683-gene cilia set**: human 8 sig (all ↓); mouse 223 (bidirectional, incidental)
- Both used the **authors' processed counts** (raw reads embargoed). Results verified **byte-reproducible**.

## Environment / gotchas
- Project `.venv`: pydeseq2, gseapy, pandas, sklearn, matplotlib, openpyxl, fastapi, uvicorn
- `~/miniforge3`: hisat2, subread(featureCounts), samtools
- Salmon binary + human index were in `/tmp` → **WIPED on WSL restart** (rebuild if needed)
- `/tmp` wipes on WSL restart → persist to `results/`; conda tool wrappers need their bin on PATH (align.py handles it); recent 2026 GEO series embargo raw reads → use processed counts

## Open threads
1. **Pan-cancer cilia markers** — TCGA via UCSC **Xena** (log2 TPM) or **recount3** (counts); decide all-cancers vs epithelial subset → per-cancer tumor-vs-normal on cilia set → recurrence heatmap
2. **`cilia_aml/`** — Cilia-in-AML; need AML-vs-normal data (TCGA-LAML has **no normals** → use GEO AML+normal or GTEx blood). Candidate seen: `SRX348003xx` (miR-551b AML/KG1a study)
3. Optional: organism-aware enrichment (mouse libraries); unify PCA methods

_To continue after `/compact` or in a fresh session: read this file + the repo._
