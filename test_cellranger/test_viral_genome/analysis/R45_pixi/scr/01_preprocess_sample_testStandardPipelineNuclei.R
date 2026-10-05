# AIM ---------------------------------------------------------------------
# Preprocessing of the samples from the cellranger count runs of test_viral_genome (count_viral/ = host+EBV reference, count_host/ = host-only reference).
# Processes each sample from both runs (map over `in_id_samples` and `runs`); pick the samples in the parameters below.
# This is the sample run on the sample nuclei dataset

# load libraries ----------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scuttle)
library(hdf5r)
library(scCustomize)

# custom functions --------------------------------------------------------


# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# parameters --------------------------------------------------------------
# cellranger runs to process: "viral" (host+EBV reference) and "host" (host-only reference)
runs <- c("viral", "host")
# runs <- c("viral")

# sample folder inside count_viral/ and count_host/, one of:
# GSM4796271_LCL_777_B958, GSM4796272_LCL_777_M81, GSM4796273_LCL_461_B958,
# GSM3596320_GM18502, GSM3596321_GM12878, GSM3596322_GM12878-GM18502_MIXTURE,
# connect_5k_pbmc_NGSC3_ch1_gex_1 (negative control)
# in_id_sample <- "connect_5k_pbmc_NGSC3_ch1_gex_1"
# in_id_sample <- "GSM4796271_LCL_777_B958"
# in_id_sample <- "GSM3596322_GM12878-GM18502_MIXTURE"
# in_id_sample <- "GSM5695152_TX1242_0_GEX"
in_id_samples <- c("GSM5695152_TX1242_0_GEX","GSM5695158_TX1242_8_GEX")

# define organism
in_id_org <- "9606"

# test_viral_genome/ folder (relative to the project root analysis/R45_pixi)
proj_dir <- "/beegfs/scratch/ric.cosr/pedrini.edoardo/pipeline/test_out/GSE189141"

# viral GTF used to build the combined reference: its genes are tagged as viral
viral_gtf <- "/beegfs/scratch/ric.absinta/ric.absinta/reference/GRCh38-2024-A_plus_EBV/ref/viral.gtf"

stopifnot(all(runs %in% c("viral", "host")))

# viral genes ---------------------------------------------------------------
# genes of the viral GTF (gene_name, as used for the matrix rownames;
# Seurat turns "_" into "-" in feature names). The host-only run has none, so it is 0 there.
viral_genes <- read_lines(viral_gtf) %>%
  str_extract('gene_name "[^"]+"') %>%
  str_remove_all('gene_name |"') %>%
  na.omit() %>%
  unique() %>%
  str_replace_all("_", "-")

walk(in_id_samples,function(in_id_sample){
  # define the output
  # dir.create("out/object", recursive = TRUE, showWarnings = FALSE)
  # dir.create("out/table", recursive = TRUE, showWarnings = FALSE)
  
  message("sample id: ", in_id_sample)
  message("id organism: ", in_id_org)
  
  # processing --------------------------------------------------------------
  # read one cellranger run of the sample, build the object with the QC metadata, save it
  # run <- "viral"
  tryCatch(walk(runs, function(run) {
    input_target <- file.path(proj_dir, run, "results/cellranger/merged",in_id_sample)
    h5_file <- file.path(input_target, "outs", "filtered_feature_bc_matrix.h5")
    out_id_object <- file.path("out/object", paste0("01_", run, "_", in_id_sample, "_obj_preQC_testStandardPipelineNuclei.rds"))
    out_id_meta <- file.path("out/table", paste0("01_", run, "_", in_id_sample, "_meta_preQC_testStandardPipelineNuclei.tsv"))
    
    message("run: ", run)
    message("Loading CellRanger output from: ", h5_file)
    message("output rds: ", out_id_object)
    message("output meta: ", out_id_meta)
    
    # Check if file exists (good practice for debugging)
    if (!file.exists(h5_file)) {
      stop(paste("The h5 file was not found at:", h5_file))
    }
    
    # load the input file ---------------------------------------------------
    # read in the matrix
    # This function can handle both default cellranger or cellbender output
    data <- Read_CellBender_h5_Mat(file_name = h5_file)
    
    # standard processing ---------------------------------------------------
    # crete the object
    datasc <- CreateSeuratObject(counts = data, project = in_id_sample,
                                 # potentially parametrize the values below
                                 # do not filter any genes before the integration
                                 min.cells = 0,
                                 min.features = 200)
    datasc$run <- run
    
    # add the mefadata accordin to the organism
    if(in_id_org == "9606"){
      # add the metadata
      datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
      datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
      datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
    } else if(in_id_org == "10090"){
      # add the metadata
      datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^mt-")
      datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
      datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^Hb[^(p)]")
    } else {
      # decide whether to run something or nothing
    }
    
    # -----------------------------------------------------------------------
    # viral load
    viral_features <- intersect(viral_genes, rownames(datasc))
    message("viral genes in the matrix: ", length(viral_features), " of ", length(viral_genes))
    
    if (length(viral_features) > 0) {
      datasc$percent.viral <- PercentageFeatureSet(datasc, features = viral_features)
      datasc$nCount_viral <- colSums(LayerData(datasc, layer = "counts")[viral_features, , drop = FALSE])
    } else {
      datasc$percent.viral <- 0
      datasc$nCount_viral <- 0
    }
    
    # -----------------------------------------------------------------------
    # add the filtering variable based on the adaptive threshold multivalue
    stats <- cbind(log10(datasc@meta.data$nCount_RNA),
                   log10(datasc@meta.data$nFeature_RNA),
                   datasc@meta.data$percent.mt)
    
    outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
    multi.outlier <- isOutlier(outlying, type = "higher")
    datasc$discard_multi <- as.vector(multi.outlier)
    
    # add the filtering variable based on the adaptive threshold single values
    high_QC_mito <- isOutlier(datasc@meta.data$percent.mt, type="high", log=TRUE)
    QC_features <- isOutlier(datasc@meta.data$nFeature_RNA, type="both", log=TRUE)
    
    datasc$discard_single <- high_QC_mito | QC_features
    
    # save output -----------------------------------------------------------
    saveRDS(datasc,out_id_object)
    write_tsv(datasc@meta.data %>% rownames_to_column("barcodes"),out_id_meta)
  }), error = function(e) message("FAILED ", in_id_sample, ": ", conditionMessage(e)))
})
