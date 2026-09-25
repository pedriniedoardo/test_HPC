#!/bin/bash
#SBATCH --job-name=cellranger_host
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=8  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="cellranger_host.err"   # submit_counts.sh overrides these with one log per sample
#SBATCH --output="cellranger_host.out"

# Step 2b - same call as 02_cellranger_count_viral.sh, against the host-only reference ($HOST_REF),
# as a baseline to compare with the host+virus run. Does not need steps 00/01.
# One sample per job. Submit all samples from test_viral_genome/:  bash scr/submit_counts.sh scr/02b_cellranger_count_host.sh
set -euo pipefail

source /beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test_viral_genome/scr/config.sh

FASTQ_DIR=${1:?"usage: sbatch $0 <sample FASTQ folder> <short id>"}
ID=${2:?"usage: sbatch $0 <sample FASTQ folder> <short id>"}
SAMPLE=$(basename "$FASTQ_DIR")   # folder is named like the FASTQ prefix
set --   # `. activate` would otherwise receive this script's arguments and fail

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_cellranger100;

cd "$PROJ"
mkdir -p "$PROJ/count_host"   # cellranger creates the run folder but not its parent

cellranger count --id=$ID \
        --transcriptome=$HOST_REF \
        --fastqs=$FASTQ_DIR \
        --sample=$SAMPLE \
        --localcores=8 \
        --create-bam=false \
        --localmem=64 \
        --output-dir=$PROJ/count_host/$ID
