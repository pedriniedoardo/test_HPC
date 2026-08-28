# AIM ---------------------------------------------------------------------
# test script to check pixi and renv integration

# renv integration --------------------------------------------------------

# to load the packages
source(".Rprofile")

# alternative implementation to be set in the rule

# params:
#   # Pass the renv library path from the config
#   renv_lib = config["renv_library_path"],
# 
# shell:
#     """
#         # Set R_LIBS to point directly to the renv library cache.
#         # This overrides R's default library search paths.
#         export R_LIBS={params.renv_lib}
#         
#         # Now run the R script
#         Rscript --vanilla scr/001_import_snakemake.R
#         """
#

# in the config specify the following
# renv_library_path: "renv/library/linux-rocky-9.5/R-4.5/x86_64-conda-linux-gnu"

# libraries ---------------------------------------------------------------
library(Seurat)
library(harmony)
library(tidyverse)

# Snakemake integation ----------------------------------------------------
# define the input
# input_id <- snakemake@input$inp
# message("input: ", input_id)

# define the output
# save the filtered dataset
output_id <- snakemake@output$out_file
message("output: ", output_id)

# wrangling ---------------------------------------------------------------
# save the ouput
writeLines(capture.output(sessionInfo()), con = output_id)

