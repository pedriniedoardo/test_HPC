#!/bin/sh
#SBATCH --job-name=test-pixi
#SBATCH --output=test_02.out
#SBATCH --error=test_02.err
#SBATCH --time=INFINITE
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --partition=workq
#SBATCH --mail-type=END
#SBATCH --mail-user=pedrini.edoardo@hsr.it

echo "job started"
date

# activate the pixi env, conda-activate style
eval "$(pixi shell-hook --manifest-path /beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_pixi_cluster/test_sbatch)"

python -c "import sys; print('python:', sys.version); print('executable:', sys.executable); print('prefix:', sys.prefix)"

date
