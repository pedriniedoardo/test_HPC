#!/bin/sh
#SBATCH --job-name=test-pixi
#SBATCH --output=test_04.out
#SBATCH --error=test_04.err
#SBATCH --time=INFINITE
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --partition=workq
#SBATCH --mail-type=END
#SBATCH --mail-user=pedrini.edoardo@hsr.it

# relies on pixi.toml being in the directory sbatch was submitted from, no --manifest-path
pixi run bash -c '
echo "job started"
date
python -c "import sys; print(\"python:\", sys.version); print(\"executable:\", sys.executable); print(\"prefix:\", sys.prefix)"
date
'
