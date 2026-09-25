#!/bin/bash
# Shared settings for 00_prepare_viral_inputs.sh, 01_mkref_host_viral.sh, 02_cellranger_count_viral.sh
# Sourced (not executed). Absolute paths on purpose: sbatch runs a copy of the script from a spool dir.

PROJ=/beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test_viral_genome

# ---- host reference (same as test100/scr/00_cellranger.sh) ----
HOST_REF=/beegfs/scratch/ric.absinta/ric.absinta/reference/refdata-gex-GRCh38-2024-A
HOST_TAG=GRCh38-2024-A

# ---- organism to add ----
# Download the FASTA + GTF yourself from UCSC GenArk (see README section 2) and put them in $PROJ/download/.
# Plain or .gz both work. They must be the same assembly version, and contig names must match
# (GTF column 1 = FASTA ">" name). EBV RefSeq = GCF_002402265.1 (contig NC_007605.1).
VIRAL_FASTA=$PROJ/download/GCF_002402265.1.fa.gz
VIRAL_GTF=$PROJ/download/GCF_002402265.1_ASM240226v1.ncbiGene.gtf.gz
VIRUS_NAME=EBV               # label used in the reference folder name only; [A-Za-z0-9_-]

# ---- outputs ----
REF_DIR=$PROJ/ref
COMBINED_NAME=${HOST_TAG}_plus_${VIRUS_NAME}
COMBINED_REF=$REF_DIR/$COMBINED_NAME      # final reference passed to cellranger count --transcriptome

# ---- samples to count (02, 02b) ----
# Dataset folders, each with one subfolder per sample named like its FASTQs: <folder>/<folder>_S1_L00*_R{1,2}_001.fastq.gz.
# submit_counts.sh submits one 02/02b job per subfolder.
SAMPLES_DIRS=(
    # $PROJ/sample_dataset/GSE158275/grouped     # 3 EBV+ lymphoblastoid cell lines
    $PROJ/sample_dataset/GSE126321/grouped     # 3 EBV+ samples
    # $PROJ/sample_dataset/connect_5k_pbmc       # 10x 5k PBMC (test100 run), negative control
)
