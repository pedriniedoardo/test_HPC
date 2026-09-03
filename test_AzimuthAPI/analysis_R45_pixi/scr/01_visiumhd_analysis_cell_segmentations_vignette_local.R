# REF ---------------------------------------------------------------------
# Analysis and visualization of Visium HD data with cell segmentations in Seurat
# Source: https://github.com/satijalab/seurat/blob/HEAD/vignettes/visiumhd_analysis_cell_segmentations.Rmd
# this is mainly to test the local panhuman annotation

# libraries ---------------------------------------------------------------
library(Seurat)
library(patchwork)
library(tidyverse)
library(hdf5r)
library(arrow)
library(sf)
library(reticulate)

# 2. Point reticulate at the local panhumanpy Python environment -------------
# Pin RETICULATE_PYTHON directly to the env's interpreter (must be set before reticulate initializes any Python session).
# This is more robust than use_condaenv(name) on this machine, since that relies on `conda env list` to discover environments, which does not reliably see micromamba envs.
# Adjust the path if your pan-human-azimuth env lives elsewhere (check with: micromamba env list).
Sys.setenv(RETICULATE_PYTHON = "./.conda-pan-human-azimuth/bin/python")

# Confirm reticulate actually loaded THIS interpreter and that it can see panhumanpy inside it. 
py_config()
# might expose token
# Sys.getenv()
Sys.getenv(c("RETICULATE_PYTHON", "CONDA_PREFIX", "CONDA_DEFAULT_ENV"))

cat("Active Python:", py_config()$python, "\n")
cat("panhumanpy version:", as.character(import("panhumanpy")$`__version__`), "\n")

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
kidney_obj_local <- ANNotate(kidney_obj,
                             output_mode = "minimal",
                             assay = "Spatial.Polygons",
                             model_version = "v1")

saveRDS(kidney_obj_local,file = "../out/object/01_kidney_obj_local.rds")
