# REF ---------------------------------------------------------------------
# Analysis and visualization of Visium HD data with cell segmentations in Seurat
# Source: https://github.com/satijalab/seurat/blob/HEAD/vignettes/visiumhd_analysis_cell_segmentations.Rmd

# ---- Overview -----------------------------------------------------------
# Traditional methods in spatial transcriptomics use bins typically ranging 2-16 um in size to approximate spatial locations of transcripts; however, binning approaches lose single-cell resolution.
# Released in June 2025, 10x Genomics' Space Ranger cell segmentation workflow supports nucleus and cell segmentation for Visium HD and Visium HD 3' H&E samples, allowing transcripts to be assigned at the cell level.
#
# The segmentation approach involves capturing gene expression at subcellular resolution (~2 um) and performing nuclei detection and expansion, which enables the delineation of individual cell boundaries. This allows much finer profiling of features, with gene counts reported on a per-cell basis rather than aggregated in 8x8um or 16x16um squares. For more information, see documentation from 10x:
# https://www.10xgenomics.com/support/software/space-ranger/latest/algorithms-overview/image-processing#cell-segmentation
#
# With Seurat v5.4+, we introduce support for these new high-resolution datasets from Space Ranger, including the ability to:
#  - Load, process, and visualize spatial gene expression with segmentation overlays,
#  - Represent cell centroids and polygonal boundaries as distinct datatypes, and
#  - Perform clustering and spatial domain identification at multiple resolutions within the same object.
#
# Together, these developments power segmentation-level spatial transcriptomics analysis in Seurat.
#
# This vignette demonstrates how to use Seurat to analyze Visium HD cell segmentations. The typical analysis pipeline for these data is identical to those introduced in previous spatial vignettes; as such, here we aim to highlight the advantages of performing analyses at the cell level rather than at the spot or bin level, especially for complex and heterogeneous tissue architectures.

# First, we load Seurat and the other packages necessary for this vignette. Core functionality in support of Visium HD cell segmentation data was initially introduced in Seurat v5.4.0 & SeuratObject v5.3.0.

# libraries ---------------------------------------------------------------
library(Seurat)
library(patchwork)
library(tidyverse)
library(hdf5r)
library(arrow)
library(sf)
library(ComplexHeatmap)
library(viridis)

# ---- Load data ------------------------------------------------------------
# In this vignette, we will use a Visium HD dataset of the human kidney (cortex) analyzed with Space Ranger 4.0.1.
# The data is available for download from the 10x data repository here:
# https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-human-kidney-ffpe-v4

# How do we set up the data directory?
# First, create a directory to hold the data.
# mkdir -p visium_hd_human_kidney
# cd visium_hd_human_kidney

# Then, go to the "Output and supplemental files" tab in the link to the repository above to access the core
# output files. The "batch download" option gives curl and wget commands to download the data via the terminal.
# wget https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Kidney_FFPE/Visium_HD_Human_Kidney_FFPE_binned_outputs.tar.gz
# wget https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Kidney_FFPE/Visium_HD_Human_Kidney_FFPE_segmented_outputs.tar.gz
# wget https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Kidney_FFPE/Visium_HD_Human_Kidney_FFPE_spatial.tar.gz
# wget https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Kidney_FFPE/Visium_HD_Human_Kidney_FFPE_barcode_mappings.parquet
# wget https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Kidney_FFPE/Visium_HD_Human_Kidney_FFPE_feature_slice.h5
# wget https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Kidney_FFPE/Visium_HD_Human_Kidney_FFPE_metrics_summary.csv
# wget https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Kidney_FFPE/Visium_HD_Human_Kidney_FFPE_molecule_info.h5

# Next, extract the binned outputs (expression matrices for 2um, 8um, and 16um resolution bins), segmented outputs (segmentation-based expression data and segmentation polygons), and spatial data (H&E tissue images and spatial coordinate mappings).
# tar -xzf Visium_HD_Human_Kidney_FFPE_binned_outputs.tar.gz
# tar -xzf Visium_HD_Human_Kidney_FFPE_segmented_outputs.tar.gz
# tar -xzf Visium_HD_Human_Kidney_FFPE_spatial.tar.gz

# After downloading the data, we load it into Seurat using the `Load10X_Spatial()` function. As introduced in previous vignettes, this function reads in the output of the Space Ranger `count` pipeline. As of Space Ranger 4.0+, the count pipeline also automatically runs the `segment` pipeline for Visium HD and Visium HD 3' data. We specify which data to load with the `bin.size` argument. By including `"polygons"`, we load the segmentation-based expression data where transcripts are assigned to individual cells based on the nuclei detection and cell boundary expansion method. In this vignette, we also load the 8um-binned data for comparison.

# More about the `Load10X_Spatial` arguments used here:
# The image.name argument was also introduced in Seurat v5.4. Currently, users specify either "tissue_lowres_image.png" (default) or "tissue_hires_image.png". The high-resolution image is 10 times larger, so it is useful for generating higher quality visualizations; however, users should also consider object size & analysis needs when deciding what to load. The documentation for `Load10X_Spatial` provides more details on other arguments that may be of interest to users who want to customize their workflows.

path_to_data <- "../data/sample_kidney_VisiumHD/"

kidney_obj <- Load10X_Spatial(data.dir = path_to_data,
                              bin.size = c(8, "polygons"),
                              image.name = "tissue_hires_image.png")

# We can inspect the images and their associated assays to verify the data stored in the object:
Images(kidney_obj)

kidney_obj[["slice1.008um"]]
kidney_obj[["slice1.polygons"]]

# ---- Bins vs. segmentations ------------------------------------------------
# One of the key advantages of the segmentation approach is the ability to analyze and visualize data at single-cell resolution. Before we plot the data on the whole tissue image, we can first define regions of interest in the tissue image for focused analysis later.

# In order to do so, we begin by inspecting the whole tissue image and its coordinate ranges.
tissue_image <- SpatialDimPlot(kidney_obj,
                               images = "slice1.008um",
                               alpha = 0,
                               image.scale = "hires") +
                    ggtitle("Tissue image") + NoLegend()

# Define layer to customize axis text
theme_layer <- theme(axis.text.x = element_text(size = 10, angle = 45),
                                    axis.text.y = element_text(size = 10),
                                    axis.title.x = element_text(size = 10),
                                    axis.title.y = element_text(size = 10))

tissue_image_with_axes <- tissue_image +
                              scale_x_continuous(n.breaks = 20) +
                              scale_y_reverse(n.breaks = 20) +
                              theme_layer

tissue_image_with_axes

# The approach we take here is to extract tissue coordinates and subset our data to the same coordinate ranges for direct comparison. (For a purely visual comparison--without subsetting or creating a new FOV--with ggplot2 v4.0.0+ users can specify the desired `xlim` and `ylim` directly in a `coord_fixed()` layer on the plot of the whole tissue image.) By using `Crop()` (https://satijalab.github.io/seurat-object/reference/Crop.html), we can define a new zoomed field of view (FOV), stored in the `images` slot.
#
# Here we define two different levels of zoom, based on coordinates from the tissue image above.
#
# For the first level of zoom, we will take a section of 800x400 pixels, and apply this to the segmented data.
zoom1_xlim <- c(2400, 3200)
zoom1_ylim <- c(3300, 3700)

kidney_obj[["zoom1.polygons"]] <- Crop(object = kidney_obj[["slice1.polygons"]],
                                              x = zoom1_xlim,
                                              y = zoom1_ylim,
                                              scale = "hires")

# For the second level of zoom, we will take a smaller section of 400x300 pixels, and apply this to both the binned and segmented data.
zoom2_xlim <- c(2750, 3150)
zoom2_ylim <- c(3350, 3650)

kidney_obj[["zoom2.polygons"]] <- Crop(object = kidney_obj[["slice1.polygons"]],
                                              x = zoom2_xlim,
                                              y = zoom2_ylim,
                                              scale = "hires")

kidney_obj[["zoom2.008um"]] <- Crop(object = kidney_obj[["slice1.008um"]],
                                           x = zoom2_xlim,
                                           y = zoom2_ylim,
                                           scale = "hires")

# To see the regions captured in these zooms, we plot just the tissue:
tissue_image_zoom1 <- SpatialDimPlot(kidney_obj,
                                     images = "zoom1.polygons",
                                     alpha = 0,
                                     image.scale = "hires") +
                                  ggtitle("Zoom level 1") + NoLegend()

tissue_image_zoom2 <- SpatialDimPlot(kidney_obj,
                                     images = "zoom2.polygons",
                                     alpha = 0,
                                     image.scale = "hires") +
                                  ggtitle("Zoom level 2") + NoLegend()

tissue_image_zoom1 + tissue_image_zoom2

# Now, we can take a closer look at how bins and cell segmentations capture spatial information by comparing how they appear when overlaid on the tissue image. Let's use the second level of zoom for this comparison, since it is more focused on a smaller region of the tissue.
#
# When generating a spatial plot of segmented data, `plot_segmentations = TRUE` displays the full boundaries of each segmentation polygon, rather than just the centroids. This allows us to see the actual shape and size of each detected cell.
sdp_8um_zoom <- SpatialDimPlot(kidney_obj,
                               images = "zoom2.008um",
                               pt.size.factor = 18,
                               stroke = 0.1,
                               alpha = 0.6,
                               image.scale = "hires",
                               cols = "turquoise4") +
                    ggtitle("8µm bins") + NoLegend()

sdp_segm_zoom <- SpatialDimPlot(kidney_obj,
                                images = "zoom2.polygons",
                                plot_segmentations = TRUE,
                                stroke = 0.1,
                                alpha = 0.6,
                                image.scale = "hires",
                                cols = "turquoise4") +
                      ggtitle("Cell segmentations") + NoLegend()

sdp_8um_zoom + sdp_segm_zoom

# ---- Cell type annotation with Pan-human Azimuth ---------------------------
# We now annotate cell types in the segmentation-level data.
#
# To perform cell type annotation, we use Pan-human Azimuth (https://satijalab.org/pan_human_azimuth/), a neural network-based classifier trained on a reference scRNA-seq dataset of 23 different types of human tissue. The classifier maps query cells to a hierarchical cell ontology and provides confidence scores for each annotation.
# The original implementation of this tool is in Python as panhumanpy (https://github.com/satijalab/panhumanpy).
#
# In R, we can access Pan-human Azimuth through the AzimuthAPI package (https://github.com/satijalab/AzimuthAPI):
# if (!requireNamespace("AzimuthAPI", quietly = TRUE)) {
#   remotes::install_github("satijalab/AzimuthAPI")
# }

library(AzimuthAPI)

# ---- Annotating cell segmentations -----------------------------------------
# `AzimuthAPI::CloudAzimuth()` runs Pan-human Azimuth on a Seurat object via a cloud-based API. The object is returned with the results of cell type annotation as added cell-level metadata, and the latent embeddings of the classifier as an added dimensional reduction. For more details, see the Pan-human Azimuth R API vignette:
# https://satijalab.org/panhuman_r_vignette
#
# We set the default assay, normalize the data with standard log-normalization, and then run `CloudAzimuth()`.
DefaultAssay(kidney_obj) <- "Spatial.Polygons"

kidney_obj <- NormalizeData(kidney_obj)

# remote implementation
kidney_obj_remote <- CloudAzimuth(kidney_obj,
                                  assay = "Spatial.Polygons",
                                  model_version = "v1")

saveRDS(kidney_obj_remote,file = "../out/object/01_kidney_obj_remote.rds")

# First, we can inspect how the output of Pan-human Azimuth is stored within the object metadata. Note that the object contains both 8um bins and cell segmentations; after running `CloudAzimuth` above on the `Spatial.Polygons` assay, the Azimuth metadata columns are populated for (only) cell segmentations.
library(knitr)
kidney_obj_remote@meta.data %>% head()
kidney_obj_remote@meta.data %>% tail()

kidney_obj_remote@meta.data %>% mutate(test = is.na(nCount_Spatial.Polygons)) %>% group_by(test) %>% summarise(n = n())

# check the levels before filtering
kidney_obj_remote@meta.data %>%
  filter(!is.na(nCount_Spatial.Polygons)) %>%
  group_by(final_level_labels) %>%
  summarise(n = n())

# We see that in addition to cell type annotations organized hierarchically, Pan-human Azimuth provides confidence scores for each annotation; lower confidence scores indicate uncertainty in hierarchy assignments. We can set a confidence threshold to retain only high-confidence annotations for downstream tasks.
kidney_obj_hc <- subset(kidney_obj_remote,
                        !is.na(final_level_confidence) & final_level_confidence > 0.7)

# check the levels after filtering
kidney_obj_hc@meta.data %>%
  filter(!is.na(nCount_Spatial.Polygons)) %>%
  group_by(final_level_labels) %>%
  summarise(n = n())

# In addition to using the confidence score as a QC metric for annotation accuracy, we also encourage users to examine differentially expressed genes for each annotated cell type with `AzimuthAPI::make_azimuth_QC_heatmaps()`.
#
# For example, we can look at a heatmap for cells annotated within the epithelial cell hierarchy:
kidney_azimuth_QC_heatmaps <- make_azimuth_QC_heatmaps(kidney_obj_hc)
kidney_azimuth_QC_heatmaps[["Epithelial cell_1"]]

# Additionally, the azimuth_embed dimensional reduction layer stores a 128-dimensional reduction of the Pan-human Azimuth embeddings. This is useful in downstream analyses. For example, we can visualize the cell type predictions returned by Pan-human Azimuth in embedding space:
kidney_obj_hc <- RunUMAP(kidney_obj_hc,
                         dims = 1:128,
                         reduction = "azimuth_embed",
                         reduction.name = "azimuth_umap")

p2 <- DimPlot(kidney_obj_hc,
              group.by = "final_level_labels",
              label.size = 2.5,
              label = TRUE,
              reduction = "azimuth_umap",
              repel = TRUE) + NoLegend()

# compare this with what I would have obtained from the regular analysis
# Clustering (standard Seurat). the object contains only the cells
kidney_obj_hc <- FindVariableFeatures(kidney_obj_hc) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20,
                reduction = "pca",
                graph.name = c("pca_nn", "pca_snn")) %>%
  FindClusters(graph.name = "pca_snn",
               resolution = seq(0.2, 1, by = 0.2)) %>%
  RunUMAP(return.model = TRUE,
          dims = 1:20,
          reduction = "pca",
          reduction.name = "pca_umap")

DimPlot(kidney_obj_hc,
        reduction = "pca_umap",
        group.by = "pca_snn_res.0.4",
        label = T,
        raster = TRUE,
        repel = T,
        label.size = 2.5) + NoLegend()

kidney_obj_hc@graphs

# ---- Comparing azimuth_embed-based and PCA-based clustering to Pan-human Azimuth labels ----
# How well does unsupervised clustering on each embedding (azimuth_embed vs. the standard PCA computed above) recover the cell type structure assigned by Pan-human Azimuth? We cluster on azimuth_embed (analogous to the PCA clustering above), then compare both clusterings to `final_level_labels` via a Jaccard similarity score computed on shared cell barcodes.
# Cluster on the azimuth_embed reduction (128 dims) at resolution 0.4, matching the resolution already used for the PCA-based clustering comparison above (pca_snn_res.0.4).
#
# NOTE: FindNeighbors()'s default graph.name is derived from DefaultAssay(), which is still "Spatial.Polygons" here -- without an explicit graph.name this would silently overwrite the "pca_nn"/"pca_snn" graphs from the PCA FindNeighbors() call above.
kidney_obj_hc <- kidney_obj_hc %>% 
  FindNeighbors(dims = 1:128,
                reduction = "azimuth_embed",
                graph.name = c("azimuth_nn", "azimuth_snn")) %>%
  FindClusters(graph.name = "azimuth_snn",
               resolution = 0.4)

kidney_obj_hc@graphs
kidney_obj_hc@meta.data %>% head()

# ---- Jaccard similarity between cluster assignments and cell type labels ----
# use the Jaccard score to measure the cross-cluster similarity per cell (how similar are the clusters from the query compared to the annotation derived from the reference)

# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# build the dataset for the correlation plot crossing the two annotations per cell, and compute the jaccard score per pair (factored into a helper since we need this twice: azimuth clusters vs. labels, and PCA clusters vs. labels)
compute_jaccard_matrix <- function(meta, id_ref_col, id_query_col) {
  df_crossing <- crossing(id_ref = unique(meta[[id_ref_col]]),
                          id_query = unique(meta[[id_query_col]]))

  df_jaccard_score <- pmap(list(id_ref = df_crossing$id_ref,
                                id_query = df_crossing$id_query), function(id_ref, id_query){

                                  a <- meta %>%
                                    filter(.data[[id_ref_col]] == id_ref) %>% pull(barcode)

                                  b <- meta %>%
                                    filter(.data[[id_query_col]] == id_query) %>% pull(barcode)

                                  jaccard_score <- jaccard(a, b)

                                  df <- data.frame("id_ref" = id_ref,
                                                   "id_query" = id_query,
                                                   "jaccard_score" = jaccard_score)
                                  return(df)
                                }) %>%
    bind_rows()

  # shape it as a matrix
  df_jaccard_score %>%
    pivot_wider(names_from = id_ref, values_from = jaccard_score) %>%
    column_to_rownames("id_query")
}

meta_hc <- kidney_obj_hc@meta.data %>% rownames_to_column("barcode")

# Heatmap A: azimuth_embed-based clusters vs. Pan-human Azimuth cell type labels
mat_jaccard_azimuth <- compute_jaccard_matrix(meta_hc, "final_level_labels", "azimuth_snn_res.0.4")

ht_jaccard_azimuth <- Heatmap(mat_jaccard_azimuth,
                              name = "Jaccard score",
                              col = viridis::viridis(option = "turbo", n = 20),
                              row_names_side = "right",
                              row_names_gp = gpar(fontsize = 8),
                              column_names_side = "bottom",
                              column_names_gp = gpar(fontsize = 8),
                              row_dend_reorder = FALSE,
                              column_dend_reorder = FALSE,
                              row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                              column_title = "final_level_labels",
                              row_title = "azimuth_embed clusters (res 0.4)",
                              show_column_names = TRUE,
                              show_row_names = TRUE)

p01 <- DimPlot(kidney_obj_hc,
               group.by = "final_level_labels",
               label.size = 2.5,
               label = TRUE,
               reduction = "azimuth_umap",
               repel = TRUE) + NoLegend()

p02 <- DimPlot(kidney_obj_hc,
               group.by = "azimuth_snn_res.0.4",
               label.size = 2.5,
               label = TRUE,
               reduction = "azimuth_umap",
               repel = TRUE) + NoLegend()

# grid.grabExpr(draw(x)) i used to make ad ggplot object a complexheatmap output
# also give more spece to the heatmap
(p01 + p02 + grid.grabExpr(draw(ht_jaccard_azimuth))) +
  plot_layout(widths = c(1,1,3))



# Heatmap B: PCA-based clusters vs. Pan-human Azimuth cell type labels
mat_jaccard_pca <- compute_jaccard_matrix(meta_hc, "final_level_labels", "pca_snn_res.0.4")

ht_jaccard_pca <- Heatmap(mat_jaccard_pca,
                          name = "Jaccard score",
                          col = viridis::viridis(option = "turbo", n = 20),
                          row_names_side = "right",
                          row_names_gp = gpar(fontsize = 8),
                          column_names_side = "bottom",
                          column_names_gp = gpar(fontsize = 8),
                          row_dend_reorder = FALSE,
                          column_dend_reorder = FALSE,
                          row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                          column_title = "final_level_labels",
                          row_title = "PCA clusters (res 0.4)",
                          show_column_names = TRUE,
                          show_row_names = TRUE)

p03 <- DimPlot(kidney_obj_hc,
               group.by = "final_level_labels",
               label.size = 2.5,
               label = TRUE,
               reduction = "pca_umap",
               repel = TRUE) + NoLegend()

p04 <- DimPlot(kidney_obj_hc,
               group.by = "pca_snn_res.0.4",
               label.size = 2.5,
               label = TRUE,
               reduction = "pca_umap",
               repel = TRUE) + NoLegend()

# grid.grabExpr(draw(x)) i used to make ad ggplot object a complexheatmap output
# also give more spece to the heatmap
(p03 + p04 + grid.grabExpr(draw(ht_jaccard_pca))) +
  plot_layout(widths = c(1,1,3))


# ---- Cell type annotations in spatial context ------------------------------
# Next, we can visualize where different cell types are located within the tissue. Mapping the cell type annotations onto the tissue image, we can observe the cell types present in various anatomically distinct regions (an introduction to the structure of the kidney cortex through an interactive view of normal tissue can be found via the Human Protein Atlas: https://www.proteinatlas.org/learn/dictionary/normal/kidney), such as glomeruli and renal tubules.
#
# Podocytes are highly specialized cells found within the glomeruli of the kidney. We can focus our visualization on only cells annotated as podocytes:
kidney_segm_pod <- subset(kidney_obj_hc, final_level_labels == "Podocyte")

Idents(kidney_segm_pod) <- "final_level_labels"

pod_annotations_tissue <- SpatialDimPlot(kidney_segm_pod,
                                         images = "slice1.polygons",
                                         plot_segmentations = TRUE,
                                         image.scale = "hires",
                                         cols = "cyan") + NoLegend()

pod_annotations <- SpatialDimPlot(kidney_segm_pod,
                                  images = "slice1.polygons",
                                  image.alpha = 0,
                                  plot_segmentations = TRUE,
                                  image.scale = "hires",
                                  cols = "turquoise4") + NoLegend()

pod_annotations_tissue +
  pod_annotations +
  plot_annotation(title = 'Cells w/ annotation \'Podocyte\'')

# We can pick out some additional cell types of interest to examine the spatial organization of major kidney cell populations.
cell_types_of_interest <- c("Capillary EC",
                            "Pericyte",
                            "Podocyte",
                            "Renal epithelial cell - distal tubules",
                            "Renal epithelial cell - Loop of Henle",
                            "Type A intercalated cell",
                            "Type B intercalated cell")

kidney_segm_top <- subset(kidney_obj_hc,
                          final_level_labels %in% cell_types_of_interest)

Idents(kidney_segm_top) <- "final_level_labels"

cell_type_colors <- c("#06d31e", "#006326", "#000000", "#453494", "#ffaf01", "#00fffa", "#1da5f9")
names(cell_type_colors) <- cell_types_of_interest

full_image_annotated <- SpatialDimPlot(kidney_segm_top,
                                       images = "slice1.polygons",
                                       plot_segmentations = TRUE,
                                       alpha = 0.7,
                                       stroke = 0.04,
                                       image.scale = "hires",
                                       cols = cell_type_colors) +
                            labs(fill = "Cell type")

full_image_annotated

# ---- High-resolution zoom into tissue regions ------------------------------
# Finally, zooming in again, we observe cell types annotated in and around some glomeruli.

# Define layers to customize the legend by putting it in three rows below the plot
legend_guide_layer <- list(theme(legend.position = "bottom",
                                 legend.text = element_text(size = 12),
                                 legend.title = element_blank()),
                           guides(fill = guide_legend(nrow = 3, byrow = TRUE)))

glom_zoom_annotated <- SpatialDimPlot(kidney_segm_top,
                                      images = "zoom1.polygons",
                                      group.by = "final_level_labels",
                                      plot_segmentations = TRUE,
                                      alpha = 0.7,
                                      stroke = 0.1,
                                      image.scale = "hires",
                                      cols = cell_type_colors) +
                            legend_guide_layer

glom_zoom_annotated

# In this section of kidney cortex, which includes several glomeruli surrounded by renal tubules, the predicted cell type annotations align well with the underlying tissue structure. Podocytes are confined to the glomeruli, while endothelial and perivascular cells cluster within capillary-rich regions as expected. Distinct tubular epithelial populations are distributed along surrounding tubules, consistent with the normal spatial organization of nephron segments.

# ---- More information --------------------------------------------------------
# For more information on spatial analysis workflows, see the other spatial vignettes
# (https://satijalab.org/seurat/articles/get_started_v5_new#spatial-analysis) and documentation
# (Seurat: https://satijalab.org/seurat/reference/, SeuratObject: https://satijalab.github.io/seurat-object/reference/index.html).

# Acknowledgements: thanks to Stephen Williams and the Computational Biology team at 10x for their helpful
# feedback and contributions to the code for loading Visium segmentation data.
