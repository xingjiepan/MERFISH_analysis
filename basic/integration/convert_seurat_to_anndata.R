# R script for converting a Seurat h5 file to an h5ad file
# Usage:
# Rscript --vanilla convert_seurat_to_anndata.R input_file

library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)
input_file = args[1]

Convert(input_file, dest='h5ad', overwrite=TRUE, verbose=TRUE)

