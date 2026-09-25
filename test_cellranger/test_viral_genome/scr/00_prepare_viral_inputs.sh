#!/bin/bash
# Step 0 - clean and check the FASTA + GTF you downloaded (VIRAL_FASTA, VIRAL_GTF in config.sh)
# into ref/viral.fa and ref/viral.gtf. Light job: run on the login node with `bash`.
set -euo pipefail

# Loads config.sh for VIRAL_FASTA, VIRAL_GTF, VIRUS_NAME, REF_DIR and HOST_REF.
source /beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test_viral_genome/scr/config.sh

mkdir -p "$REF_DIR"
FA=$REF_DIR/viral.fa
GTF=$REF_DIR/viral.gtf


# Checks that VIRUS_NAME contains only letters, digits, _ and -. It becomes part of the reference folder name (GRCh38-2024-A_plus_EBV), and mkref rejects other characters.
[[ "$VIRUS_NAME" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "ERROR: VIRUS_NAME must match [A-Za-z0-9_-]+" >&2; exit 1; }

# ---- 1. read the input files ------------------------------------------------------
# Stops if either file is missing or empty.
for f in "$VIRAL_FASTA" "$VIRAL_GTF"; do
    [[ -s "$f" ]] || { echo "ERROR: $f not found or empty; download it and set VIRAL_FASTA/VIRAL_GTF in config.sh" >&2; exit 1; }
done

# Decompresses them if they're .gz, or copies them as they are if not, into temporary viral.raw.fa and viral.raw.gtf.
zcat -f "$VIRAL_FASTA" > "$REF_DIR/viral.raw.fa"
zcat -f "$VIRAL_GTF"   > "$REF_DIR/viral.raw.gtf"

# Writes the two input paths to viral.source.txt. Step 01 compares this with config.sh, so it refuses to run on files prepared from an older config.
printf '%s\n' "$VIRAL_FASTA" "$VIRAL_GTF" > "$REF_DIR/viral.source.txt"

# ---- 2. clean ---------------------------------------------------------------
# FASTA → ref/viral.fa:
# It shortens each header to its first word, so >NC_007605.1 Human gammaherpesvirus 4, complete genome becomes >NC_007605.1. The GTF refers to the sequence by that one word, and Cell Ranger uses it as the chromosome name.
# It converts the sequence to uppercase (UCSC writes repeats in lowercase), removes stray spaces and Windows line endings, and rewraps lines at 60 characters.
awk '
    function flush() { if (seq != "") { for (i = 1; i <= length(seq); i += 60) print substr(seq, i, 60); seq = "" } }
    /^>/ { flush(); print $1; next }
    { gsub(/[ \t\r]/, ""); seq = seq toupper($0) }
    END { flush() }
' "$REF_DIR/viral.raw.fa" > "$FA"

# GTF → ref/viral.gtf:
# It drops comment lines and malformed lines.
# It replaces / with _ in the attributes. Some EBV genes are named like BBLF2/BBLF3, and a / in a gene name can break file paths or tools later on.
awk -F'\t' -v OFS='\t' '!/^#/ && NF >= 9 { gsub("/", "_", $9); print }' "$REF_DIR/viral.raw.gtf" > "$GTF"

# It deletes the temporary raw files.
rm "$REF_DIR/viral.raw.fa" "$REF_DIR/viral.raw.gtf"

# It writes viral.lengths.tsv with each sequence's name and length, NC_007605.1  171823, which the checks below use.
awk '/^>/ { if (id) print id "\t" len; id = substr($1, 2); len = 0; next } { len += length($0) } END { print id "\t" len }' "$FA" > "$REF_DIR/viral.lengths.tsv"

# ---- 3. sanity checks ---------------------------------------------------------
# Sanity checks
# Each check stops the script with an error message rather than letting mkref fail hours later:

# 1. Name clash with human: if a viral sequence had the same name as a human chromosome, the combined reference would be broken.
if cut -f1 "$HOST_REF/fasta/genome.fa.fai" | grep -qxFf <(cut -f1 "$REF_DIR/viral.lengths.tsv"); then
    echo "ERROR: a viral contig name already exists in the host genome" >&2; exit 1
fi

# 2. The GTF names a sequence that isn't in the FASTA: the usual sign that the two files come from different sources or assembly versions.
awk -F'\t' '
    NR == FNR { len[$1] = $2; next }
    !($1 in len)   { print "ERROR: GTF contig " $1 " not in viral FASTA (contig names differ?)" > "/dev/stderr"; bad = 1 }
    ($5 > len[$1]) { print "ERROR: GTF feature beyond contig end: " $1 ":" $4 "-" $5 > "/dev/stderr"; bad = 1 }
    $3 == "exon" { nexon++ }
    END { if (!nexon) { print "ERROR: no exon feature in the GTF; mkref counts only exons" > "/dev/stderr"; bad = 1 } exit bad }
' "$REF_DIR/viral.lengths.tsv" "$GTF"

# 3. A gene runs past the end of the sequence: another sign of mismatched versions.
echo "OK: $FA"; cat "$REF_DIR/viral.lengths.tsv"

# 4. No exon lines at all: mkref builds genes only from exons, so nothing would be counted. This would catch NCBI-style CDS-only GTFs if they had no exons at all.
echo "OK: $GTF"
awk -F'\t' '$3 == "exon" { n++; match($9, /gene_id "[^"]+"/); g[substr($9, RSTART + 9, RLENGTH - 10)]; s[$7]++ }
            END { k = 0; for (x in g) k++; print k " genes, " n " exons (+ strand: " s["+"] + 0 ", - strand: " s["-"] + 0 ")" }' "$GTF"
