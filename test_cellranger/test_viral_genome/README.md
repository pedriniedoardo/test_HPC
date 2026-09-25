# Adding a virus (EBV) to the Cell Ranger reference

Cell Ranger only counts genes that are in its reference. The human reference has no viral genes, so viral reads are lost. The fix: add the virus's genome sequence and gene annotation to the human reference, then count again. Viral genes then show up as extra rows in the count matrix.

## Step by step

Run everything from `test_viral_genome/`. Each step names what it does and why it is needed. The sections below give the details.

**Step 1. Choose the genome of the organism to add.**
Open [NCBI Genomes](https://www.ncbi.nlm.nih.gov/datasets/genome/), search the organism (e.g. "EBV"), and copy the accession of the row tagged **NCBI RefSeq** (`GCF_...`), see [section 1](#1-find-the-genome-accession-on-ncbi). **The accession pins one exact genome version. Its sequence and gene annotation then match base for base. RefSeq is NCBI's curated reference, one per organism.**

**Step 2. Download the FASTA and GTF.**
Build the UCSC GenArk folder URL from the accession, open it, and download the FASTA and the GTF into `download/`, see [section 2](#2-download-the-fasta-and-gtf-from-ucsc-genark).
**A reference needs two files. The FASTA is the DNA sequence the reads are aligned to. The GTF says where each gene lies on that sequence, which lets Cell Ranger turn an aligned read into a count for a gene. Both must come from the same assembly, so the contig names and coordinates match.**

**Step 3. Point the settings at the files and prepare them.**
Edit `scr/config.sh`: set `VIRAL_FASTA`, `VIRAL_GTF` (paths to the downloaded files, `.gz` is fine) and `VIRUS_NAME=EBV`. Then:
```bash
bash scr/00_prepare_viral_inputs.sh      # seconds, fine on the login node
```
**All three scripts read their paths and names from `config.sh`, so a new organism only changes these lines. `VIRUS_NAME` is only used to name the output folder. Step 00 copies the files into `ref/viral.fa` and `ref/viral.gtf`, cleans the names, and checks that the two files agree (same contig names, genes inside the sequence, `exon` lines present). It prints the number of genes found; for EBV it is 86.**

**Step 4. Build the combined human + virus reference.**
```bash
sbatch scr/01_mkref_host_viral.sh        # hours: re-indexes the whole human genome
```
**`cellranger count` accepts only one reference. The script appends the viral FASTA and GTF to the human ones, so the virus becomes one more chromosome. `cellranger mkref` then builds the STAR alignment index. This is the step described in the 10x tutorial, which appends GFP to zebrafish the same way. Output: `ref/GRCh38-2024-A_plus_EBV/`. Check it worked: `grep -x NC_007605.1 ref/GRCh38-2024-A_plus_EBV/star/chrName.txt` prints the viral contig.**

**Step 5. Count with the new reference (test run on a sample dataset).**
This step is only a test of the new reference, not part of building it. It runs on a public sample dataset: three 10x 3' scRNA-seq samples of EBV-positive lymphoblastoid cell lines from GEO [GSE158275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158275), in `sample_dataset/GSE158275/grouped/`. A negative control sits in `sample_dataset/connect_5k_pbmc/`: links to the 10x `connect_5k_pbmc_NGSC3_ch1_gex_1` PBMC FASTQs (the `test100` run). Healthy blood has almost no EBV-infected cells, so viral counts there should be near zero. For your own data, list your dataset folders in `SAMPLES_DIRS` in `config.sh`.
```bash
bash scr/submit_counts.sh scr/02_cellranger_count_viral.sh    # host+EBV reference -> count_viral/<id>/
bash scr/submit_counts.sh scr/02b_cellranger_count_host.sh    # optional baseline: host only -> count_host/<id>/
```
**This is a normal `cellranger count` with `--transcriptome` set to the combined reference. Reads from viral RNA now align to the viral contig and are counted. Viral genes appear as extra rows (features) in the count matrix, next to the human genes. `submit_counts.sh` submits one independent job per sample folder in the `SAMPLES_DIRS` folders (`config.sh`: the three GSE158275 LCL samples and the PBMC control), with logs in `logs/`. One sample by hand:**
`sbatch scr/02_cellranger_count_viral.sh <path to sample FASTQ folder> <short id>`.
**Check it worked: `zcat count_viral/GSM4796271_LCL_777_B958/outs/filtered_feature_bc_matrix/features.tsv.gz | grep BZLF1`.**

**Step 6. Read the result.**
Load the matrix as usual (Seurat/scanpy) and look at the viral genes' counts per cell. Cells with viral UMIs are candidate infected cells. Read the EBV caveats below before interpreting zeros: some viral RNAs are not captured or not counted.

**For another virus**, repeat steps 1–5 with its files and a new `VIRUS_NAME`. Human is not re-downloaded, but step 4 must be re-run because the index changes. For where to get the files, see [where to find the FASTA and GTF](#another-virus-or-organism-where-to-find-the-fasta-and-gtf).

## 1. Find the genome accession on NCBI

1. Open [NCBI Genomes](https://www.ncbi.nlm.nih.gov/datasets/genome/) and search the organism, e.g. `EBV`
   or `Human gammaherpesvirus 4`.
2. The result table has one row per genome assembly. "EBV" returns ~586 rows, almost all single-strain submissions (`GCA_...`). Look for the row tagged **NCBI RefSeq**: the reference genome curated by NCBI, with gene annotation. For EBV: `ASM240226v1`, **`GCF_002402265.1`**.
3. Click it. The genome page (`https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002402265.1/`) confirms the accession and shows the assembly name (`ASM240226v1`) and the sequence it contains (just oen chromosome since it is a virus `NC_007605.1`, 171,823 bp, B95-8 strain with the missing part filled in from Raji).

The accession has three parts:
```
GCF _ 002402265 . 1
 │      │         └─ version (changes if the sequence is corrected)
 │      └─ 9-digit number
 └─ GCF = RefSeq (NCBI curated)   GCA = GenBank (as submitted)
```

## 2. Download the FASTA and GTF from UCSC GenArk

Search the accession on the [UCSC assembly search](https://genome.ucsc.edu/assemblySearch.html) to confirm GenArk has it. It opens the genome browser, not the files.

[UCSC GenArk](https://hgdownload.soe.ucsc.edu/hubs/) mirrors NCBI assemblies and stores each one in a folder whose path is built from the accession. We use it rather than NCBI because its GTF has `exon` lines for every gene. NCBI's own `*_genomic.gtf.gz` has them for only 16 of 94 EBV genes, and `mkref` counts only exons. The sequence is the same as NCBI's.

**a. Build the folder URL from the accession.** Split the 9-digit number into three groups of three:
```
accession        GCF_002402265.1
                  │    │  │  │
folder path      GCF/002/402/265/GCF_002402265.1/

URL = https://hgdownload.soe.ucsc.edu/hubs/  +  GCF/002/402/265/GCF_002402265.1/
    = https://hgdownload.soe.ucsc.edu/hubs/GCF/002/402/265/GCF_002402265.1/
```

**b. Open the folder in a browser.** The two files you need:
```
GCF_002402265.1/
├── GCF_002402265.1.fa.gz                          <- FASTA (the sequence)
└── genes/
    └── GCF_002402265.1_ASM240226v1.ncbiGene.gtf.gz  <- GTF (the gene annotation)
```
- The FASTA is `<accession>.fa.gz`, directly in the folder.
- The GTF is in `genes/`. Take the file ending in `.gtf.gz`. If there are several, prefer `ncbiGene` or `ncbiRefSeq` (NCBI annotation) over `augustus` (a computer prediction). The file name adds the assembly name (`_ASM240226v1`), so copy it from the listing rather than guessing it.
- Ignore the rest (`.2bit`, `.bb`, `bbi/`, `ixIxx/`, `hub.txt`, ...): those are genome-browser files.
- If the folder does not exist, GenArk does not have the assembly. If `genes/` is missing, it has no annotation. See [another virus](#another-virus-or-organism-where-to-find-the-fasta-and-gtf) for NCBI instead.

**c. Download them into `download/`** (right-click → copy link address in the listing, then `wget`):
```bash
mkdir -p download && cd download
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/002/402/265/GCF_002402265.1/GCF_002402265.1.fa.gz
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/002/402/265/GCF_002402265.1/genes/GCF_002402265.1_ASM240226v1.ncbiGene.gtf.gz
cd ..
```
The same folder URL from the terminal, for any accession:
```bash
ACC=GCF_002402265.1
N=${ACC#*_}; N=${N%.*}                                  # 002402265
echo https://hgdownload.soe.ucsc.edu/hubs/${ACC%%_*}/${N:0:3}/${N:3:3}/${N:6:3}/$ACC/
```

**d. Set them in `scr/config.sh`** (these are the current values):
```bash
VIRAL_FASTA=$PROJ/download/GCF_002402265.1.fa.gz
VIRAL_GTF=$PROJ/download/GCF_002402265.1_ASM240226v1.ncbiGene.gtf.gz
VIRUS_NAME=EBV        # only used in the name of the combined reference
```

**Alternative route: NCBI FTP.** On the NCBI genome page
(`https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002402265.1/`), click **FTP**. It opens the assembly folder, whose name adds the assembly name to the accession:
```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/
```
Download the two `_genomic` files: `GCF_002402265.1_ASM240226v1_genomic.fna.gz` (FASTA) and `GCF_002402265.1_ASM240226v1_genomic.gtf.gz` (GTF). Ignore the rest: `.gff`/`.gbff` are other formats,
`cds_from_genomic`/`rna_from_genomic`/`protein` are per-gene sequences, not the genome.
The FASTA is the same sequence as UCSC's. The GTF is the same annotation, but in NCBI's file only 16 of 96 EBV genes have `exon` lines (the rest are `CDS` only), and `mkref` counts only exons. So for EBV the UCSC GTF gives more genes (86). This can differ for other organisms: count the genes with exons before choosing (command in [another virus](#another-virus-or-organism-where-to-find-the-fasta-and-gtf)).

`bash scr/00_prepare_viral_inputs.sh` then writes into `ref/`:

- `viral.fa`: the sequence, one contig (`NC_007605.1`). UCSC writes repeats in lowercase; step 00
  converts to uppercase.
- `viral.gtf`: the gene locations, one line per exon: contig, start, end, strand, and attributes
  (`gene_id`, `gene_name`). For EBV: 86 genes, 232 exons. `mkref` only uses the `exon` lines.

## 3. Build the reference and count

```bash
cd /beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_cellranger/test_viral_genome
bash   scr/00_prepare_viral_inputs.sh   # download/ files -> ref/viral.fa + ref/viral.gtf (seconds, login node)
sbatch scr/01_mkref_host_viral.sh       # human + virus -> ref/GRCh38-2024-A_plus_EBV (slow: human genome index)
bash scr/submit_counts.sh scr/02_cellranger_count_viral.sh  # one job per sample, combined reference -> count_viral/<id>/
bash scr/submit_counts.sh scr/02b_cellranger_count_host.sh  # same samples, human-only reference -> count_host/<id>/
```

`config.sh` holds every path and setting. `02` is the `test100` script with `--transcriptome` pointed
at the combined reference (and the cloud cell-annotation flags removed), run once per sample in
`SAMPLES_DIRS` by `submit_counts.sh`. `02b` is the same with the human-only reference.

## Another virus or organism: where to find the FASTA and GTF

The FASTA and the GTF must come from the **same assembly version** and use the **same contig names** (column 1 of the GTF = the `>` name in the FASTA). The GTF needs `exon` lines with `gene_id` and `transcript_id`, because `mkref` only reads exons. GFF3 files are not accepted.

**A. A virus or other small genome added to human (this pipeline).**
1. Search the organism on [NCBI Genomes](https://www.ncbi.nlm.nih.gov/datasets/genome/) and copy the accession of the **NCBI RefSeq** row (`GCF_...`). `GCA_...` rows are unreviewed submissions and often have no annotation.
2. Download the FASTA and GTF for that accession into `download/`, set `VIRAL_FASTA`, `VIRAL_GTF` and `VIRUS_NAME` in `scr/config.sh`, then run steps 00 → 01 → 02. The folder paths come from the accession (`GCF_002402265.1` → `GCF/002/402/265/`):
   - **UCSC GenArk** (first choice, as in [section 2](#2-download-the-fasta-and-gtf-from-ucsc-genark)):
     `https://hgdownload.soe.ucsc.edu/hubs/GCF/002/402/265/GCF_002402265.1/`
     → `GCF_002402265.1.fa.gz` + `genes/*.gtf.gz`. Not every accession is there, and some have no `genes/`.
   - **NCBI FTP** (if GenArk has no GTF; the FTP button on the NCBI genome page): `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/`
     → folder `GCF_002402265.1_ASM240226v1/` → `*_genomic.fna.gz` + `*_genomic.gtf.gz`
     (or `datasets download genome accession GCF_002402265.1 --include genome,gtf`).
     Before using NCBI's GTF, count the genes that have exon lines:
     `zcat file.gtf.gz | awk -F'\t' '$3=="exon"' | grep -o 'gene_id "[^"]*"' | sort -u | wc -l`.
     For EBV it has exons for only 16 of 94 genes.
     If many genes only have `CDS` lines, turn the CDS lines into exons:
     `awk -F'\t' -v OFS='\t' '$3=="CDS"{$3="exon"; $8="."; print}'`.
   - **A single sequence** (e.g. a strain with no assembly, or a transgene like GFP/Cre): take the FASTA from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) (*Send to → File → FASTA*), shorten the header to one word, and write the GTF by hand: one `exon` line per gene/ORF. Use coordinates from the GenBank record, or 1..length for a single transcript, as in the 10x GFP example:
     `echo -e 'GFP\tunknown\texon\t1\t922\t.\t+\t.\tgene_id "GFP"; transcript_id "GFP"; gene_name "GFP";'`
3. Check the files in IGV (load the FASTA as genome, the GTF as track). 00 also checks that contigs match and that coordinates stay inside the contig.

**B. A whole new host organism (e.g. mouse, macaque, zebrafish), as in the 10x tutorial.** Use
[Ensembl](https://www.ensembl.org/info/data/ftp/index.html): `fasta/<species>/dna/*.dna.primary_assembly.fa.gz`
(or `*.dna.toplevel.fa.gz` if there is no primary assembly) and `gtf/<species>/<Species>.<assembly>.<release>.gtf.gz`
from the **same release**. Ensembl GTFs have `gene_biotype`, so filter them first:
`cellranger mkgtf in.gtf out.gtf --attribute=gene_biotype:protein_coding --attribute=gene_biotype:lncRNA ...`.
Then point `HOST_REF` at the result of `cellranger mkref`, or at a
[10x pre-built reference](https://www.10xgenomics.com/support/software/cell-ranger/downloads) (human, mouse).
The pre-built references are already filtered, which is why this pipeline skips `mkgtf`.

## Caveats for EBV

- Many EBV genes overlap. Reads falling in exons of two genes on the same strand are dropped: 10.8% of `+`-strand exon bases (2.9% on `-`) lie under two or more genes. Counting whole gene spans including introns gives 67%, mostly the long spliced EBNA transcripts, but only exons matter here.
- Few 3' UTRs are annotated (17 transcripts). Reads primed on a polyA tail can fall outside every gene. Most EBV exons stop at the stop codon, but 10x 3' reads come mostly from the 3' UTR after it, so they align just past the annotated gene and are not counted. Expect EBV counts to be underestimates; comparing one gene across cells is fine, comparing levels between EBV genes is not.
  *Possible workaround (not applied, the GTF is used as is):* extend the last exon of every transcript by ~300–500 bp in the 3' direction (to the right on `+`, to the left on `-`, not past the contig ends), between steps 00 and 01. Tested on EBV with 500 bp: `mkref` accepts it and keeps all 86 genes. The cost is more overlap between neighbouring genes, whose shared reads are dropped (exon bases under ≥2 genes: `+` 10.8% → 21.2%, `-` 2.9% → 15.8%). Whether it gains more reads than it loses is untested; check on EBV-positive data before adopting it.
- EBER RNAs are not polyadenylated, so 3' capture will under-represent them.
- LMP-2A/2B are not in the UCSC GTF, so they have no row in the matrix. Their transcripts span the point where the circular genome is cut into the linear sequence (exon 1 near the end, exons 2–9 at 58–1,682), which the UCSC conversion drops. A missing LMP-2 row does not mean LMP-2 is not expressed.

## Check the result

```bash
grep -x NC_007605.1 ref/GRCh38-2024-A_plus_EBV/star/chrName.txt                        # viral contig is in the index
zcat count_viral/GSM4796271_LCL_777_B958/outs/filtered_feature_bc_matrix/features.tsv.gz | grep BZLF1   # viral gene is a row
```

The GSE158275 samples are EBV-transformed LCLs, so viral genes should have counts in most cells.
`count_host/` is the same data without EBV in the reference (no viral rows by construction).

## Optional: data with EBV

Lymphoblastoid cell lines carry EBV in every cell. The eLife paper below deposited 10x scRNA-seq of
them (GEO GSE158275, plus GSE126321) and used the same approach (`NC_007605` as an extra chromosome).
To see what the custom reference adds, run `cellranger count` twice on the same FASTQs, once with
`refdata-gex-GRCh38-2024-A` and once with the combined reference, and compare the EBV counts (the
first has none by construction). Not checked: whether raw reads are deposited (the GEO pages could not
be opened from here).

## Status

- Step 00 was run for EBV (with the old automatic download; the local-file version was tested on the same two files and gives identical `ref/viral.fa` and `ref/viral.gtf`). Steps 01 and 02 have **not** been run on the full human+EBV reference.
- `mkref` was tested on EBV alone (seconds): all 86 genes reach the STAR index. With `--output-dir`,
  cellranger 10 writes the reference **into that directory**, even though its final message prints
  `./<genome>`. `01` therefore uses `--output-dir=$COMBINED_REF`.
- Viral genes show gene type `MissingGeneType`, because the GTF has no `gene_biotype`. This does
  not affect counting.

## References

- 10x Genomics, [Build a Custom Reference for Cell Ranger (mkref)](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr)
- [NCBI Genomes](https://www.ncbi.nlm.nih.gov/datasets/genome/)
- eLife, [scRNA-seq of lymphoblastoid cell lines](https://elifesciences.org/articles/62586)
