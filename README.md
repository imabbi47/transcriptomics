# transcriptomics — ingest stage

Fetch RNA-seq data from public archives by **accession** and turn it into
analysis-ready FASTQ + a samplesheet that downstream steps consume.

This is **Stage 0 (ingest)** of a larger, organism-agnostic RNA-seq pipeline
(`fetch → QC/trim → quantify → differential expression → enrichment → report`).
It is deliberately small and self-contained.

## What it does

```
accession  ──►  resolve to runs  ──►  download FASTQ  ──►  samplesheet.csv + runs.json
(SRR/SRX/SRP,    (ENA + NCBI)         (ENA direct, or         (input for the
 ERR/DRR, PRJ,                         SRA Toolkit)            next stage)
 GSE/GSM)
```

- **Accepts** any of: SRA/ENA/DDBJ runs (`SRR…/ERR…/DRR…`), experiments (`SRX…`),
  samples (`SRS…`), studies (`SRP…`), BioProjects (`PRJNA…`), and **GEO**
  series/samples (`GSE…/GSM…`).
- **Resolves** metadata via the **ENA portal API** (direct FASTQ URLs) and **NCBI
  E-utilities** (the reliable path from GEO → SRA runs).
- **Downloads** from **ENA** when possible (already-gzipped FASTQ, MD5-verified),
  falling back to the **SRA Toolkit** (`prefetch` + `fasterq-dump`) otherwise.
- **Emits** `samplesheet.csv` (sample, run, fastq_1, fastq_2, strandedness,
  organism, layout) and `runs.json` (full provenance).

## Install

```bash
conda env create -f environment.yml
conda activate transcriptomics
pip install -e .
```

The Python code is **stdlib-only**; the conda env just provides the external
tools the downloader uses (`sra-tools`, `pigz`, `curl`).

## Usage

```bash
# Inspect what an accession contains — metadata only, no download:
transcriptomics resolve GSE60450
transcriptomics resolve SRP043510

# Download to ./data and write a samplesheet:
transcriptomics fetch SRR1039508 -o data

# A whole GEO series, paired or single auto-detected:
transcriptomics fetch GSE60450 -o data

# Quick test: only the first 100k reads via the SRA route:
transcriptomics fetch SRR1039508 --method sra --max-spots 100000 -o data

# Dry run: resolve + write planned samplesheet/runs.json, download nothing:
transcriptomics fetch GSE60450 --dry-run -o data
```

Also runnable without installing: `python -m transcriptomics resolve SRR1039508`.

### Options (`fetch`)

| Flag | Meaning |
|------|---------|
| `-o, --outdir` | output directory (default `./data`) |
| `--method {auto,ena,sra}` | download route (default `auto`) |
| `--threads N` | threads for `fasterq-dump` (default: auto = cores − 2) |
| `--max-spots N` | download only the first N reads (quick tests; SRA route) |
| `--overwrite` | re-download even if files exist |
| `--keep-intermediates` | keep `.sra`/work files |
| `--dry-run` | resolve + metadata only |

### Environment variables

- `NCBI_API_KEY` — raises NCBI rate limit from 3 → 10 requests/s (recommended).
- `NCBI_EMAIL` — contact address (NCBI E-utilities etiquette).

## Output layout

```
data/
├── SRR1039508_R1.fastq.gz
├── SRR1039508_R2.fastq.gz
├── runs.json          # full resolved metadata
└── samplesheet.csv    # input contract for the next stage
```

## Tests

```bash
python -m pytest          # or: python -m unittest discover tests
```

## Status / roadmap

- [x] Stage 0 — ingest: `fetch` / `resolve` (stdlib core)
- [ ] Stage 1 — QC + trimming (`fastp`)  ← next
- [x] Stage 2 — quantification: `quantify` (Salmon)
- [x] Stage 3 — differential expression: `de` (pyDESeq2)
- [x] Stage 4 — enrichment: `enrich` (gseapy) + `report` (single-file HTML)
- [x] helper — `geo-design` (parse a GEO series matrix into a design table)

## Downstream analysis

Needs the heavier stack (`de` / `enrich`); `quantify` / `report` / `geo-design` are stdlib:

```bash
pip install -e '.[analysis]'

# quantify FASTQ against a Salmon index
transcriptomics quantify --salmon salmon --index idx --r1 R1.fq.gz --r2 R2.fq.gz --outdir quant

# differential expression from a count matrix + design sheet
transcriptomics de --counts counts.tsv.gz --design design.csv \
    --factor treatment --ref Control --alt Treated --covariate "batch,sex" --outdir de

# enrichment, then one self-contained HTML report
transcriptomics enrich --results de/de_results.csv --outdir de --label "Treated vs Control"
transcriptomics report --outdir de --title "My study" --out report.html
```

Validated end-to-end on **GSE334363** (mouse; Sildenafil vs Control → cholesterol/immune signature)
and **GSE336399** (human; Cancer vs Benign → immune suppression) — see `analysis/`.
