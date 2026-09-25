#!/bin/bash
#SBATCH --job-name=mkref_viral
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=8  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="mkref_viral.err"
#SBATCH --output="mkref_viral.out"

# Step 1 - host + virus FASTA/GTF -> cellranger mkref. Submit from test_viral_genome/:  sbatch scr/01_mkref_host_viral.sh
set -euo pipefail

source /beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test_viral_genome/scr/config.sh

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_cellranger100;

[[ -s "$REF_DIR/viral.fa" && -s "$REF_DIR/viral.gtf" ]] || { echo "Run scr/00_prepare_viral_inputs.sh first" >&2; exit 1; }
[[ "$(cat "$REF_DIR/viral.source.txt" 2>/dev/null)" == "$(printf '%s\n' "$VIRAL_FASTA" "$VIRAL_GTF")" ]] \
    || { echo "ref/viral.fa|gtf were not made from the VIRAL_FASTA/VIRAL_GTF in config.sh (stale?); re-run scr/00_prepare_viral_inputs.sh" >&2; exit 1; }
[[ ! -e "$COMBINED_REF" ]] || { echo "$COMBINED_REF already exists; remove it or change VIRUS_NAME" >&2; exit 1; }

# `awk 1` re-emits every line with a trailing newline, so a host file lacking one cannot glue onto the viral record
awk 1 "$HOST_REF/fasta/genome.fa" "$REF_DIR/viral.fa"  > "$REF_DIR/host_plus_viral.fa"
awk 1 "$HOST_REF/genes/genes.gtf" "$REF_DIR/viral.gtf" > "$REF_DIR/host_plus_viral.gtf"

# With --output-dir, mkref (cellranger 10) writes the reference itself (fasta/ genes/ star/ reference.json) into that
# directory, despite its final message printing ./<genome> (checked with an EBV-only test build).
cd "$REF_DIR"
cellranger mkref \
        --genome=$COMBINED_NAME \
        --fasta=$REF_DIR/host_plus_viral.fa \
        --genes=$REF_DIR/host_plus_viral.gtf \
        --nthreads=8 \
        --memgb=60 \
        --output-dir=$COMBINED_REF

[[ -s "$COMBINED_REF/reference.json" ]] || { echo "mkref finished but $COMBINED_REF/reference.json is missing" >&2; exit 1; }
grep -qx "$(cut -f1 "$REF_DIR/viral.lengths.tsv" | head -1)" "$COMBINED_REF/star/chrName.txt" \
    || { echo "viral contig missing from $COMBINED_REF/star/chrName.txt" >&2; exit 1; }
echo "OK: viral contig in STAR index"
# the concatenated inputs are copied into the reference (fasta/genome.fa, genes/genes.gtf.gz); drop the 4.6 GB intermediates
rm "$REF_DIR/host_plus_viral.fa" "$REF_DIR/host_plus_viral.gtf"
