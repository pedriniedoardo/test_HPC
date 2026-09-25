#!/bin/bash
# Submits one independent sbatch job per sample folder in $SAMPLES_DIRS (config.sh), running the given count script.
# Run from test_viral_genome/ on the login node:
#   bash scr/submit_counts.sh scr/02_cellranger_count_viral.sh     # host+virus reference -> count_viral/<id>/
#   bash scr/submit_counts.sh scr/02b_cellranger_count_host.sh     # host-only reference  -> count_host/<id>/
# Logs: logs/<script>_<id>.out|err
set -euo pipefail

SCRIPT=${1:?"usage: bash $0 <count script, e.g. scr/02_cellranger_count_viral.sh>"}
[[ -f "$SCRIPT" ]] || { echo "$SCRIPT not found" >&2; exit 1; }

source /beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test_viral_genome/scr/config.sh

mkdir -p "$PROJ/logs"
for dataset in "${SAMPLES_DIRS[@]}"; do
for dir in "$dataset"/*/; do
    dir=${dir%/}
    SAMPLE=$(basename "$dir")
    # GSM4796271__LCL_777_B958__Homo_sapiens__... -> GSM4796271_LCL_777_B958 (the full name is too long for --id);
    # names without "__" (e.g. connect_5k_pbmc_NGSC3_ch1_gex_1) are used as they are
    ID=$SAMPLE
    if [[ $SAMPLE == *__* ]]; then
        REST=${SAMPLE#*__}
        ID=${SAMPLE%%__*}_${REST%%__*}
    fi
    NAME=$(basename "$SCRIPT" .sh)_$ID
    echo -n "$ID: "
    
    # Mock run: prints the command to the screen but does NOT execute sbatch
    echo sbatch --job-name=\"$NAME\" --output=\"$PROJ/logs/$NAME.out\" --error=\"$PROJ/logs/$NAME.err\" \"$SCRIPT\" \"$dir\" \"$ID\"
done
done
