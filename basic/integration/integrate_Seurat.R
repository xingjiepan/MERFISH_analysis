# R script for integrating a query Seurat dataset to a reference Seurat dataset.
# Usage:
# Rscript --vanilla integrate_Seurat.R ref_file query_file output_path cell_type_col drop_gene
# Arguments:
#   ref_file: A Seurat H5 file for the reference dataset.
#   query_file: A Seurat H5 file for the query dataset.
#   output_path: Path to save the output files.
#   cell_type_col: The column name of cell type labels. 
#   drop_gene: The gene to dropout during integration.


library(Seurat)
library(patchwork)
library(SeuratDisk)

# Define the parameters
args = commandArgs(trailingOnly=TRUE)
ref_file = args[1]
query_file = args[2]
output_path = args[3]
cell_type_col = args[4]

if (length(args) >=5){
    drop_gene = args[5]
} else {
    drop_gene = ''
}

# Load input data
dR <- LoadH5Seurat(ref_file)
dR <- AddMetaData(dR, factor(c('reference')), col.name = 'source')
dQ <- LoadH5Seurat(query_file)
dQ <- AddMetaData(dQ, factor(c('query')), col.name = 'source')

d_list <- list(dR, dQ)

# normalize and identify variable features for each dataset independently
d_list <- lapply(X = d_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 600)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = d_list)
features = features[features != drop_gene]
anchors_t <- FindTransferAnchors(reference = dR, query = dQ, reduction='cca', features=features)

predictions_counts <- TransferData(anchorset = anchors_t, weight.reduction = 'cca', refdata=attr(attr(dR, 'assay')[[1]], 'counts'))
imputed_Q <- CreateSeuratObject(counts = predictions_counts)

predictions_class_label <- TransferData(anchorset = anchors_t, weight.reduction = 'cca', refdata=dR@meta.data[[cell_type_col]])
imputed_Q <- AddMetaData(imputed_Q, predictions_class_label[1], col.name = 'predicted_cell_type')
imputed_Q <- AddMetaData(imputed_Q, predictions_class_label['prediction.score.max'], col.name = 'predicted_cell_type_proba')

SaveH5Seurat(imputed_Q, filename = paste(output_path, 'imputation.h5Seurat', sep='/'), overwrite = TRUE)
Convert(paste(output_path, 'imputation.h5Seurat', sep='/'), dest = "h5ad", overwrite = TRUE)
