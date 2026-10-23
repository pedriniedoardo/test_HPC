#!/bin/bash
# Submit step 01 or 02 for the 3 new Flex samples (one independent job per sample).
# Usage (login node, from test_flex/):
#   bash scr/submit_all.sh 01     # bamtofastq for all samples
#   bash scr/submit_all.sh 02     # cellranger multi for all samples (needs config/multi_config_<ID>.csv with the right probe barcodes)
#   bash scr/submit_all.sh 01 18a_24hr_LPC   # only the listed sample(s)
set -euo pipefail

STEP=${1:?"usage: bash $0 <01|02> [sample ...]"}
shift
IDS=("$@")
(( ${#IDS[@]} )) || IDS=(18a_24hr_EtOH 18a_24hr_LPC 18a_8week_LPC 18a_8week_EtOH)

for ID in "${IDS[@]}"; do
    case $STEP in
        01)
            [[ -s data/${ID}_sample_alignments.bam ]] || { echo "data/${ID}_sample_alignments.bam missing, skipped" >&2; continue; }
            sbatch --export=ALL,ID=$ID -J bamtofastq_$ID -o bamtofastq_$ID.out -e bamtofastq_$ID.err scr/01_bamtofastq_flex.sh ;;
        02)
            [[ -f config/multi_config_$ID.csv ]] || { echo "config/multi_config_$ID.csv missing, skipped" >&2; continue; }
            sbatch -J cellranger_$ID -o cellranger_$ID.out -e cellranger_$ID.err scr/02_cellranger_multi_flex.sh config/multi_config_$ID.csv $ID ;;
        *) echo "step must be 01 or 02" >&2; exit 1 ;;
    esac
done
