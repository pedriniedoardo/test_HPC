# test_flex

Test of re-running Cell Ranger on a 10x Flex dataset (4 samples, see table below) that is only available
on SRA as Cell Ranger BAMs. The BAM is converted back to FASTQs with 10x `bamtofastq`, then `cellranger multi` is run.

```
data/       input BAM (+ README with the download command)
config/     cellranger multi config csv
scr/        config.sh + sbatch scripts + submit_all.sh (submits 01/02 for the 3 new samples)
analysis/   outputs (fastq/, bamtofastq/, cellranger_multi/)
```

## Source of the data

Public Flex dataset on SRA. Each file is the submitted Cell Ranger BAM (original format), hosted on the SRA S3 bucket:

| Sample (`ID`) | Run | BAM URL |
|---|---|---|
| `18a_8week_EtOH` (default) | [SRR34678188](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR34678188&display=data-access) | https://sra-pub-src-2.s3.amazonaws.com/SRR34678188/18a_8week_EtOH_sample_alignments.bam.1 |
| `18a_24hr_EtOH` | SRR34678190 | https://sra-pub-src-2.s3.amazonaws.com/SRR34678190/18a_24hr_EtOH_sample_alignments.bam.1 |
| `18a_24hr_LPC` | SRR34678189 | https://sra-pub-src-1.s3.amazonaws.com/SRR34678189/18a_24hr_LPC_sample_alignments.bam.1 |
| `18a_8week_LPC` | SRR34678187 | https://sra-pub-src-2.s3.amazonaws.com/SRR34678187/18a_8week_LPC_sample_alignments.bam.1 |

```bash
# sample SRA file from Flex technology (example: first sample)
wget -c -O 18a_8week_EtOH_sample_alignments.bam \
   "https://sra-pub-src-2.s3.amazonaws.com/SRR34678188/18a_8week_EtOH_sample_alignments.bam.1"
```

(all four commands are in `data/README`; the `.1` suffix of the remote name is dropped with `-O`. The local name must be
`<ID>_sample_alignments.bam`, which is what `scr/config.sh` expects.)

## Running the other samples

The scripts are per sample: `scr/config.sh` takes the sample name from the `ID` variable (default `18a_8week_EtOH`), and
derives the BAM, FASTQ and bamtofastq paths from it. For each of the 3 new samples (`<ID>` = `18a_24hr_EtOH`, `18a_24hr_LPC`, `18a_8week_LPC`):

1. download the BAM (command in `data/README`; the 3 BAMs are already in `data/`);
2. `bash scr/submit_all.sh 01` submits `01_bamtofastq_flex.sh` for the 3 samples (one job each, logs `bamtofastq_<ID>.{out,err}`;
   samples whose BAM is missing are skipped; add a sample name to run only that one, e.g. `bash scr/submit_all.sh 01 18a_24hr_LPC`);
3. copy `config/multi_config_18a_8week_EtOH.csv` to `config/multi_config_<ID>.csv` and replace the sample name in `fastq_id`, the
   `fastqs` path and `sample_id`; **set `probe_barcode_ids` for that sample** (see step 3: do not reuse `BC011|BC012|BC013`).
   This needs the FASTQs from step 2, so 01 and 02 are not chained;
4. `bash scr/submit_all.sh 02` submits `02_cellranger_multi_flex.sh` for the 3 samples (logs `cellranger_<ID>.{out,err}`;
   samples without a config csv are skipped). Equivalent single call:
   `sbatch -J cellranger_<ID> -o cellranger_<ID>.out -e cellranger_<ID>.err scr/02_cellranger_multi_flex.sh config/multi_config_<ID>.csv <ID>`.

Outputs go to `analysis/{bamtofastq,fastq}/<ID>/` and `analysis/cellranger_multi/<ID>/`. The 3 new samples are not run yet.

## What is in the BAM

- Per-sample BAM from `cellranger multi` (Cell Ranger 7.2.0, probe set `Chromium Human Transcriptome Probe Set v1.0.1`).
  It holds the reads of this sample only, not the whole multiplexed library, and only reads that Cell Ranger kept
  (reads without a probe match are not in it). Regenerated FASTQs are therefore a subset of the original data and
  counts will be close to, but not identical with, the original run.
- The header has the `10x_bam_to_fastq` `@CO` lines, so `bamtofastq` works.

## Steps (run everything from `test_flex/`)

### 0. Download the BAM (login node / screen, not in sbatch)

Command in `data/README`:

```bash
cd data
wget -c -O 18a_8week_EtOH_sample_alignments.bam \
  "https://sra-pub-src-2.s3.amazonaws.com/SRR34678188/18a_8week_EtOH_sample_alignments.bam.1"
```

Wait until the download has finished: step 1 stops if `samtools quickcheck` finds the BAM truncated.

### 1. Create the pixi env (once)

In `test_flex/` (where `pixi.toml` must live):

```bash
pixi init -c conda-forge -c bioconda   # if there is no pixi.toml yet; bioconda is required
pixi add samtools 10x_bamtofastq
```

Both tools are on bioconda only: with conda-forge alone `pixi add samtools` fails
(`No candidates were found for samtools *`). If `pixi.toml` already exists, add the channel first:
`pixi project channel add bioconda`. Check that both tools work:
`pixi run bamtofastq --help`, `pixi run samtools --version`.
(`bamtofastq` v1.4.1 also ships inside the cellranger install, in `lib/bin/`.)

### 2. BAM -> FASTQ

```bash
sbatch scr/01_bamtofastq_flex.sh
```

Settings are in `scr/config.sh`. The script activates the pixi env with `pixi shell-hook`, checks the BAM
(`samtools quickcheck`, `@CO` header), runs `bamtofastq`, and symlinks the result into one flat folder:

- `analysis/bamtofastq/18a_8week_EtOH/<read group>/...` raw `bamtofastq` output
- `analysis/fastq/18a_8week_EtOH/18a_8week_EtOH_S1_L00N_R{1,2}_001.fastq.gz` input for step 3 (symlinks; lanes L003 and L004 as in the BAM)

Check the `ls -l` at the end of `bamtofastq_flex.out` (written in the folder you submit from).

### 3. Edit the multi config

`config/multi_config_18a_8week_EtOH.csv` holds the reference (`refdata-gex-GRCh38-2020-A`), the probe set
(v1.0.1, GRCh38-2020-A, copied from the cellranger 10.0.0 install to `/beegfs/scratch/ric.absinta/ric.absinta/reference/probe_sets/`), the FASTQ folder and the sample's probe barcodes.

**Probe barcodes for `18a_8week_EtOH`: `BC011|BC012|BC013`.** The other 3 samples must be determined the same way (they may use different barcodes). They are not in the BAM header; they were read from R2 of the regenerated FASTQs:
the 8 bases at position 68 (0-based) of R2 were matched against Cell Ranger's own table
(`cellranger-10.0.0/lib/python/cellranger/barcodes/translation/probe-barcodes-fixed-rna-profiling-rna.txt`).
On 3 subsets (200k-500k reads each, from lanes 3 and 4) the split was BC011 ~36-37%, BC013 ~32%, BC012 ~29-30%, ~2.5% no match,
so this sample uses three probe barcodes. Several barcodes for one sample are separated by `|` in `probe_barcode_ids`
(syntax from memory of the 10x multi config docs: not verified on Cell Ranger 10: check it against the docs before the first run).
`fastq_id` must stay equal to the FASTQ prefix (`18a_8week_EtOH`).

Optional: GEO has the original filtered/raw matrices for this sample (GSM9125400, `*_sample_filtered_feature_bc_matrix.h5`).
The number of barcodes in the filtered file is the original cell number (usable as `force-cells,N` in `[gene-expression]`),
or run without `force-cells` and subset the raw matrix to the original barcodes. Not done yet.

To use the newer 2024-A reference and probe set (v1.1.0), change both `reference` and `probe-set` together.

### 4. cellranger multi

```bash
sbatch scr/02_cellranger_multi_flex.sh config/multi_config_18a_8week_EtOH.csv 18a_8week_EtOH
```

`cellranger multi` has no `--transcriptome` flag: the reference is read from the csv.

## Outputs

In `analysis/cellranger_multi/<ID>/outs/` (shown for `18a_8week_EtOH`):

- `per_sample_outs/18a_8week_EtOH/count/sample_filtered_feature_bc_matrix/` (+ `.h5`): counts per cell (cells x probe genes), the file to load in Seurat/Scanpy
- `multi/count/raw_feature_bc_matrix/`: all barcodes, before splitting by sample
- `per_sample_outs/18a_8week_EtOH/web_summary.html`, `metrics_summary.csv`: QC

Exact folder names are the standard `multi` ones: verify after the run.
