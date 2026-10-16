#!/bin/bash
# Shared settings for 01_bamtofastq_flex.sh and 02_cellranger_multi_flex.sh. Sourced (not executed).
# Absolute paths on purpose: sbatch runs a copy of the script from a spool dir.

PROJ=/beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test_flex

# pixi env (pixi.toml in $PROJ) providing bamtofastq + samtools, activated by 01 with `pixi shell-hook`
PIXI=/home/pedrini.edoardo/.pixi/bin/pixi

# ---- input BAM (SRR34678188, see data/README) ----
BAM=$PROJ/data/18a_8week_EtOH_sample_alignments.bam
ID=18a_8week_EtOH            # run id / FASTQ prefix / sample_id (no spaces)

# ---- outputs ----
FASTQ_DIR=$PROJ/analysis/fastq/$ID          # flat folder of <ID>_S1_L00N_R{1,2}_001.fastq.gz, input of cellranger multi
BTF_DIR=$PROJ/analysis/bamtofastq/$ID       # raw bamtofastq output (one subfolder per read group)
MULTI_DIR=$PROJ/analysis/cellranger_multi
