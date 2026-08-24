#!/usr/bin/env bash
# End-to-end heavy stage for Cilia+AML, as ONE detached job:
#   ensure gencode -> tx2gene -> (salmon index || aria2 download) -> quant x8 -> gene matrix
set +e
PROJ=/home/abhi/github-projects/brainstroming_transcriptomics
cd "$PROJ" || exit 1
export PATH="$HOME/miniforge3/bin:$PATH"
source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null
conda activate salmon 2>/dev/null || source activate salmon 2>/dev/null
SRABIN=$(ls -d "$PROJ"/cilia_aml/data/sratoolkit.*/bin 2>/dev/null | head -1)
export PATH="$SRABIN:$PATH"

D=cilia_aml/data
REF=$D/ref
IDX=$REF/salmon_index
LOG=$D/pipeline.log
: > "$LOG"
say(){ echo "[$(date +%H:%M:%S)] $*" >> "$LOG"; }
RUNS="SRR40166669 SRR40166670 SRR40166671 SRR40166672 SRR40166673 SRR40166674 SRR40166675 SRR40166676"
BASEURL="https://sra-pub-run-odp.s3.amazonaws.com/sra"

mkdir -p "$D/fastq" "$D/sra" "$D/tmp" "$D/quant" "$REF"
say "start; salmon=$(salmon --version 2>&1 | tr -d '\n'); aria2=$(aria2c --version 2>/dev/null | head -1)"

# GENCODE transcriptome (resume past flaky ~30 MB drops until gzip verifies whole)
GFA=$REF/gencode.transcripts.fa.gz
GURL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
for i in $(seq 1 60); do
  if gzip -t "$GFA" 2>/dev/null; then say "gencode complete ($(wc -c < "$GFA") bytes)"; break; fi
  say "gencode resume $i (have $(wc -c < "$GFA" 2>/dev/null || echo 0) bytes)"
  curl -sL -C - -o "$GFA" "$GURL" --max-time 300
  sleep 3
done
gzip -t "$GFA" 2>/dev/null || say "FATAL: gencode transcriptome incomplete"

if [ ! -s "$REF/tx2gene.tsv" ]; then
  say "building tx2gene"
  zcat "$GFA" | grep "^>" | sed 's/^>//' | awk -F'|' '{print $1"\t"$6}' > "$REF/tx2gene.tsv"
  say "tx2gene lines=$(wc -l < "$REF/tx2gene.tsv")"
fi

build_index(){
  if [ -f "$IDX/info.json" ]; then say "A: index exists, skip"; return; fi
  say "A: building salmon index (k=31)"
  salmon index -t "$GFA" -i "$IDX" -k 31 -p 4 --gencode >> "$LOG" 2>&1
  say "A: index build rc=$?"
}

download_reads(){
  command -v pigz >/dev/null && local PIGZ=pigz || local PIGZ=gzip
  for r in $RUNS; do
    if ls "$D"/fastq/${r}*.fastq.gz >/dev/null 2>&1; then say "B: $r present, skip"; continue; fi
    local sra="$D/sra/$r.sra"
    local done=0
    for a in 1 2 3 4 5 6 7 8; do
      aria2c -x16 -s16 --max-connection-per-server=16 --max-tries=0 --retry-wait=3 \
        --timeout=60 --file-allocation=none --summary-interval=0 --console-log-level=warn \
        -c --dir="$D/sra" -o "$r.sra" "$BASEURL/$r/$r" >> "$LOG" 2>&1
      # aria2 deletes the .aria2 control file only once the transfer is complete
      if [ -s "$sra" ] && [ ! -f "$sra.aria2" ]; then done=1; break; fi
      say "B: $r partial after try $a ($(du -h "$sra" 2>/dev/null | cut -f1)); resuming"
      sleep 3
    done
    if [ "$done" != 1 ]; then say "B: $r DOWNLOAD FAILED"; rm -f "$sra" "$sra.aria2"; continue; fi
    say "B: $r downloaded ($(ls -lh "$sra" | awk '{print $5}')); fasterq-dump"
    if fasterq-dump -e 4 -t "$D/tmp" -O "$D/fastq" "$sra" >> "$LOG" 2>&1; then
      $PIGZ -f "$D"/fastq/${r}*.fastq >> "$LOG" 2>&1
      say "B: $r done ($(ls -lh "$D"/fastq/${r}*.fastq.gz 2>/dev/null | awk '{print $5}'))"
    else
      say "B: $r FASTERQ FAILED rc=$?"
    fi
    rm -f "$sra" "$sra.aria2"
  done
  say "B: all reads done ($(ls "$D"/fastq/*.fastq.gz 2>/dev/null | wc -l)/8)"
}

build_index & PIDA=$!
download_reads & PIDB=$!
wait $PIDA
wait $PIDB
say "index + downloads complete"

for r in $RUNS; do
  fq=$(ls "$D"/fastq/${r}*.fastq.gz 2>/dev/null | head -1)
  if [ -z "$fq" ]; then say "C: $r fastq missing, skip"; continue; fi
  if [ -f "$D/quant/$r/quant.sf" ]; then say "C: $r quant exists, skip"; continue; fi
  say "C: salmon quant $r"
  salmon quant -i "$IDX" -l A -r "$fq" -p 4 -o "$D/quant/$r" >> "$LOG" 2>&1
  mr=$(grep -i "Mapping rate" "$D/quant/$r/logs/salmon_quant.log" 2>/dev/null | tail -1)
  say "C: $r quant rc=$? $mr"
done

say "D: aggregate -> counts_gene.tsv"
"$PROJ"/.venv/bin/python cilia_aml/scripts/aggregate.py "$D/quant" "$REF/tx2gene.tsv" "$D/counts_gene.tsv" >> "$LOG" 2>&1
say "D: aggregate rc=$?"
say "PIPELINE DONE ($(ls "$D"/quant/*/quant.sf 2>/dev/null | wc -l)/8 quantified)"
