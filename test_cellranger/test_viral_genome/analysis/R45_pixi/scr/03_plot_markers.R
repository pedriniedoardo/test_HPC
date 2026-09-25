# AIM ---------------------------------------------------------------------
# Plot sets of host and viral (EBV) markers for one sample, as UMAP and as dotplot per cluster,
# for the objects of the viral run (host+EBV reference) and of the host run (host-only reference).
# Reads the objects written by scr/02_QC_sample.R.

# load libraries ----------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# parameters --------------------------------------------------------------
# sample to plot (same id as in scr/01 and scr/02)
# id_sample <- "connect_5k_pbmc_NGSC3_ch1_gex_1"
id_sample <- "GSM4796271_LCL_777_B958"

# runs to plot: "viral" (host+EBV reference) and/or "host" (host-only reference)
runs <- c("viral", "host")

# cluster column used for the dotplot (FindClusters default resolution in scr/02)
group_by <- "seurat_clusters"

# marker sets. EBV names as in the matrix: Seurat turns "_" into "-" (e.g. BSLF2_BMLF1 -> BSLF2-BMLF1).
# Genes missing from an object (e.g. EBV genes in the host run) are skipped.
list_markers_host <- list(
  B_CELL = c("MS4A1", "CD19", "CD79A", "CD79B", "PAX5"),
  LCL_ACTIVATION = c("CD40", "CD80", "CD86", "ICAM1", "EBI3", "CCL22", "TRAF1", "NFKBIA"),
  PLASMA = c("XBP1", "PRDM1", "JCHAIN", "SDC1", "IGHG1"),
  PROLIFERATION = c("MKI67", "TOP2A"),
  T_NK = c("CD3E", "CD8A", "IL7R", "NKG7", "GNLY"),
  MYELOID = c("CD14", "LYZ", "FCGR3A")
)

list_markers_viral <- list(
  EBV_LATENT = c("EBNA-1", "EBNA-2", "EBNA-3A", "EBNA-3B-EBNA-3C", "EBNA-LP", "LMP-1", "RPMS1"),
  EBV_LYTIC_IE = c("BZLF1", "BRLF1"),
  EBV_LYTIC_E = c("BMRF1", "BALF2", "BSLF2-BMLF1", "BALF5", "BHRF1"),
  EBV_LYTIC_L = c("BLLF1", "BcLF1", "BLRF2", "BALF4")
)

list_markers <- c(list_markers_host, list_markers_viral)

# per-cell summary of the viral load (from scr/01)
meta_viral <- c("percent.viral", "nCount_viral")

# dir.create("out/plot", recursive = TRUE, showWarnings = FALSE)

# plots -------------------------------------------------------------------
# for each run: read the object, then save the UMAPs and the dotplot of the markers
# run <- "viral"
walk(runs, function(run) {
  rds <- file.path("out/object", paste0("02_", run, "_", id_sample, "_obj_postQC.rds"))
  if (!file.exists(rds)) {
    stop(paste("Input not found:", rds, "- run scr/01 and scr/02 with run", run))
  }
  message("run: ", run, "  input: ", rds)
  scobj <- readRDS(rds)
  title <- paste(id_sample, run, sep = " | ")

  # keep only the markers present in this object, and drop the empty sets
  list_markers_run <- map(list_markers, ~ intersect(.x, rownames(scobj))) %>%
    keep(~ length(.x) > 0)
  missing <- setdiff(unlist(list_markers), unlist(list_markers_run))
  if (length(missing) > 0) message("  markers not in the object: ", paste(missing, collapse = ", "))

  # UMAP by cluster and viral load
  list_plot_viral <- map(meta_viral, function(x) {
    FeaturePlot(scobj, features = x, order = T, reduction = "umap", raster = F) +
      scale_color_viridis_c(option = "turbo")
  })
  p_umap_viral <- wrap_plots(c(list(DimPlot(scobj, reduction = "umap", group.by = group_by, label = T, raster = F)),
                               list_plot_viral), nrow = 1) +
    plot_annotation(title = title)
  ggsave(plot = p_umap_viral,
         file.path("out/plot", paste0("03_", run, "_", id_sample, "_UMAP_clusters_viralLoad.pdf")),
         width = 15, height = 4.5)

  # UMAP of the markers
  list_plot <- map(unlist(list_markers_run, use.names = F), function(x) {
    FeaturePlot(scobj, features = x, order = T, reduction = "umap", raster = F, pt.size = 0.2) +
      scale_color_viridis_c(option = "turbo")
  })
  n_col <- 6
  p_umap_markers <- wrap_plots(list_plot, ncol = n_col) +
    plot_annotation(title = title)
  ggsave(plot = p_umap_markers,
         file.path("out/plot", paste0("03_", run, "_", id_sample, "_UMAP_markers.pdf")),
         width = 3.5 * n_col, height = 3 * ceiling(length(list_plot) / n_col), limitsize = F)

  # dotplot of the markers per cluster
  p_dot <- DotPlot(scobj,
                   features = list_markers_run,
                   dot.scale = 6,
                   cluster.idents = T,
                   group.by = group_by) +
    RotatedAxis() +
    labs(title = title) +
    theme(strip.text = element_text(angle = 90))
  ggsave(plot = p_dot,
         file.path("out/plot", paste0("03_", run, "_", id_sample, "_DotPlot_markers.pdf")),
         width = 4 + 0.3 * length(unlist(list_markers_run)), height = 6, limitsize = F)
})
