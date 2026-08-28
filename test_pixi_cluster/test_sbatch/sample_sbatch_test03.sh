#!/bin/sh
#SBATCH --job-name=test-pixi
#SBATCH --output=test_03.out
#SBATCH --error=test_03.err
#SBATCH --time=INFINITE
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --partition=workq
#SBATCH --mail-type=END
#SBATCH --mail-user=pedrini.edoardo@hsr.it

# activate the pixi env and run the tool on the env
pixi run --manifest-path /beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_pixi_cluster/test_sbatch \
    python -c "import sys; print('python:', sys.version); print('executable:', sys.executable); print('prefix:', sys.prefix)"
