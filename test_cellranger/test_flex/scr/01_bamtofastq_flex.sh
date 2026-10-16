#!/bin/bash
#SBATCH --job-name=bamtofastq_flex
#SBATCH --account pedrini.edoardo
#SBATCH --mem=32GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=8  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="bamtofastq_flex.err"
#SBATCH --output="bamtofastq_flex.out"

# Step 1 - rebuild the FASTQs of the Flex BAM (SRR34678188) with 10x bamtofastq.
# Usage:  sbatch scr/01_bamtofastq_flex.sh      (settings in config.sh)
# Output: $FASTQ_DIR/<ID>_S1_L00N_R{1,2}_001.fastq.gz  (flat symlinks to the bamtofastq output, one lane number per read group)
set -euo pipefail

source /beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test_flex/scr/config.sh

[[ -f "$PROJ/pixi.toml" ]] || { echo "$PROJ/pixi.toml not found: create the pixi env with bamtofastq + samtools" >&2; exit 1; }
cd "$PROJ"
eval "$("$PIXI" shell-hook)"   # `pixi shell` is interactive and does not work in batch scripts

[[ -f "$BAM" ]] || { echo "$BAM not found" >&2; exit 1; }

# a truncated BAM (interrupted download) lacks the BGZF EOF marker
samtools quickcheck -v "$BAM" || { echo "$BAM is truncated/corrupt: incomplete download? see data/README" >&2; exit 1; }
samtools view -H "$BAM" | grep -q '^@CO.*10x_bam_to_fastq' || { echo "no 10x_bam_to_fastq @CO header: bamtofastq cannot be used" >&2; exit 1; }

# bamtofastq refuses to write into an existing folder
rm -rf "$BTF_DIR" "$FASTQ_DIR"
mkdir -p "$(dirname "$BTF_DIR")" "$FASTQ_DIR"

bamtofastq --nthreads=8 "$BAM" "$BTF_DIR"

# Flatten: cellranger multi wants one folder whose files start with the fastq_id.
# bamtofastq keeps the real lane in the file names (L003, L004...), so only the prefix is renamed;
# `ln -s` fails (set -e) if two read-group folders ever contain the same lane.
for f in "$BTF_DIR"/*/*.fastq.gz; do
    base=$(basename "$f")                       # bamtofastq_S1_L003_R1_001.fastq.gz
    ln -s "$f" "$FASTQ_DIR/${base/bamtofastq_/${ID}_}"
done

echo "FASTQs in $FASTQ_DIR:"
ls -l "$FASTQ_DIR"
