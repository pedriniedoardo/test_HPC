# AIM ---------------------------------------------------------------------
# Quick QC of one sample: fixed-threshold filters + doublet detection, then reprocessing.
# Reads the objects written by scr/01_preprocess_sample.R for the same sample, for both runs (map over `runs`).
# The filters below start fully open (nothing removed): tighten them after looking at the data.

# load libraries ----------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scDblFinder)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# parameters --------------------------------------------------------------
# cellranger runs to process: "viral" (host+EBV reference) and "host" (host-only reference)
runs <- c("viral", "host")
# runs <- c("viral")

# same sample as in scr/01_preprocess_sample.R
# id_sample <- "connect_5k_pbmc_NGSC3_ch1_gex_1"
id_sample <- "GSM4796271_LCL_777_B958"
# id_sample <- "GSM3596322_GM12878-GM18502_MIXTURE"

# fixed-threshold filters (open by default: nothing is removed), used for both runs
nCount_RNA_low <- 0
nCount_RNA_high <- Inf

nFeature_RNA_low <- 0
nFeature_RNA_high <- Inf

percent.ribo_low <- 0
percent.ribo_high <- 100

percent.mt_low <- 0
percent.mt_high <- 100

percent.globin_low <- 0
percent.globin_high <- 100

# remove the cells called doublets by scDblFinder
remove_doublets <- FALSE

# minimum nCount_RNA to keep a cell for the doublet identification and downstream
min_nCount_doublets <- 500

stopifnot(all(runs %in% c("viral", "host")))
message("sample id: ", id_sample)

# build the label
label <- paste("nCount-",nCount_RNA_low,"_",nCount_RNA_high,"|",
               "nFeat-",nFeature_RNA_low,"_",nFeature_RNA_high,"|",
               "pRibo-",percent.ribo_low,"_",percent.ribo_high,"|",
               "pMito-",percent.mt_low,"_",percent.mt_high,"|",
               "pHemo-",percent.globin_low,"_",percent.globin_high,"|",
               "wDbl-",remove_doublets,"|",
               "V5",sep = "")

# processing --------------------------------------------------------------
# filter one run of the sample, call doublets, reprocess and save it
walk(runs, function(run) {
  rds <- file.path("out/object", paste0("01_", run, "_", id_sample, "_obj_preQC_testStandardPipeline.rds"))
  out_id_object <- file.path("out/object", paste0("02_", run, "_", id_sample, "_obj_postQC_testStandardPipeline.rds"))
  out_id_meta <- file.path("out/table", paste0("02_", run, "_", id_sample, "_meta_postQC_testStandardPipeline.rds"))
  
  message("run: ", run)
  message("input rds: ", rds)
  message("output rds: ", out_id_object)
  message("output meta: ", out_id_meta)
  
  if (!file.exists(rds)) {
    stop(paste("Input not found:", rds, "- run scr/01_preprocess_sample.R first"))
  }
  
  # standard processing ---------------------------------------------------
  scobj <- readRDS(rds)
  
  # add the filtering label
  scobj$label <- label
  
  # add the filtering variable based on the fixed threshold
  scobj$discard_threshold <- scobj@meta.data %>%
    mutate(test = nCount_RNA < nCount_RNA_low | nCount_RNA > nCount_RNA_high |
             nFeature_RNA < nFeature_RNA_low | nFeature_RNA > nFeature_RNA_high |
             percent.ribo < percent.ribo_low | percent.ribo > percent.ribo_high |
             percent.mt < percent.mt_low | percent.mt > percent.mt_high |
             percent.globin < percent.globin_low | percent.globin > percent.globin_high) %>%
    pull(test)
  
  # preprocess the dataset before the doublet identification as recommended in:
  # https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
  # 1.5.11
  # according to the documentation the doublet identification should be run before the fine filtering
  # remove low coverage cells and proprocess to generata clusters needed for the doublet identification
  scobj <- subset(scobj,subset = nCount_RNA > min_nCount_doublets) %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters() %>%
    # do I run the UMAP ? I do not need it for the doublet identification, but can be useful in case someone wants to explore an individual sample
    RunUMAP(dims = 1:30)
  
  # run scDblFinder after filtering the low coverage cells
  sce_scobj <- scDblFinder(GetAssayData(scobj, layer="counts"), clusters=Idents(scobj))
  
  # port the resulting scores back to the Seurat object:
  scobj$scDblFinder.score <- sce_scobj$scDblFinder.score
  scobj$scDblFinder.class <- sce_scobj$scDblFinder.class
  
  # cross-tabulate doublet calls against the threshold-based discard flag
  print(table(scobj$scDblFinder.class,scobj$discard_threshold))
  
  # perform the filtering based on the fixed threshold defined, and on doublets if requested
  if(remove_doublets){
    # perform the filtering based on the fixed threshold defined and for the doublets
    scobj_filter <- subset(scobj, subset = discard_threshold == 0 & scDblFinder.class == "singlet")
  }else{
    # perform the filtering based on the fixed threshold defined only
    scobj_filter <- subset(scobj, subset = discard_threshold == 0)
  }
  message("cells kept: ", ncol(scobj_filter), " of ", ncol(scobj))
  
  # preprocess data after filtering
  scobj_filter <- scobj_filter %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:30)
  
  # save output -----------------------------------------------------------
  saveRDS(scobj_filter,out_id_object)
  write_tsv(scobj_filter@meta.data %>% rownames_to_column("barcodes"),out_id_meta)
})
