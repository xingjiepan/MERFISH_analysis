# R script for converting an h5ad file to a Seurat h5 file
# Usage:
# Rscript --vanilla convert_anndata_to_seurat.R input_file

library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)
input_file = args[1]

Convert(input_file, dest='h5seurat', overwrite=TRUE, verbose=TRUE)

