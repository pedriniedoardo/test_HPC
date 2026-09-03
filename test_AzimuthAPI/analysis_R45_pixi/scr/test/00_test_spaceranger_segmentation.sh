#!/bin/bash
#SBATCH --job-name=space
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=8  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="spaceranger_segment.err"
#SBATCH --output="spaceranger_segment.out"

echo "my job strart now" > spaceranger_segment.log;

date >> spaceranger_segment.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_spaceranger401;

spaceranger segment --id=test_kidney_sample \
    --tissue-image=Visium_HD_Human_Kidney_FFPE_tissue_image.btf \
    --localcores=8 \
    --localmem=64

date >> spaceranger_segment.log;