#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=8  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="cellranger_test_pbmc_token.err"
#SBATCH --output="cellranger_test_pbmc_token.out"

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_cellranger100;

cellranger count --id=test_pbmc_token \
        --transcriptome=/beegfs/scratch/ric.absinta/ric.absinta/reference/refdata-gex-GRCh38-2024-A \
        --fastqs=/beegfs/scratch/ric.absinta/ric.absinta/sample_data/sample_pbmcs/connect_5k_pbmc_NGSC3_ch1_fastqs/connect_5k_pbmc_NGSC3_ch1_gex_1 \
        --sample=connect_5k_pbmc_NGSC3_ch1_gex_1 \
        --localcores=8 \
        --create-bam=false \
        --tenx-cloud-token-path=/beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test100/token.txt \
        --cell-annotation-model=auto \
        --localmem=64 \
        --output-dir=/beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test100/test_pbmc_token