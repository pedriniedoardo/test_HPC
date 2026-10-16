#!/bin/bash
#SBATCH --job-name=cellranger_flex
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=8  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="cellranger_flex.err"
#SBATCH --output="cellranger_flex.out"

# Step 2 - cellranger multi (Flex) on the FASTQs made by 01_bamtofastq_flex.sh.
# Usage:  sbatch scr/02_cellranger_multi_flex.sh <multi config csv> <run id>
#   e.g.  sbatch scr/02_cellranger_multi_flex.sh config/multi_config_18a_8week_EtOH.csv 18a_8week_EtOH
# The config csv (reference, probe set, FASTQ folder, probe barcode) is written by hand: see config/multi_config_18a_8week_EtOH.csv
set -euo pipefail

source /beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test_flex/scr/config.sh

CSV=$(readlink -f "${1:?"usage: sbatch $0 <multi config csv> <run id>"}")   # absolute: sbatch runs a copy of this script
RUN_ID=${2:?"usage: sbatch $0 <multi config csv> <run id>"}
set --   # `. activate` would otherwise receive this script's arguments and fail

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_cellranger101;

[[ -f "$CSV" ]] || { echo "$CSV not found" >&2; exit 1; }
mkdir -p "$MULTI_DIR"
cd "$MULTI_DIR"

cellranger multi --id=$RUN_ID \
        --csv=$CSV \
        --localcores=8 \
        --localmem=64 \
        --output-dir=$MULTI_DIR/$RUN_ID
