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

# Only keep the genes in the query dataset for normalization
all_query_genes <- rownames(dQ)
dR_query_genes <- subset(dR, features = all_query_genes)

d_list <- list(dR_query_genes, dQ)

# Normalize and identify variable features for each dataset independently
d_list <- lapply(X = d_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst")
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = d_list)
features = features[features != drop_gene]
print('Integration using the following genes:')
print(features)

## IMPUTATION AND LABEL TRANSFER
# Find the anchors for imputation and label transfer
anchors_t <- FindTransferAnchors(reference = dR_query_genes, query = dQ, reduction='cca', features=features)

# Predict the class labels and same the predictions as a csv file 
predictions_class_label <- TransferData(anchorset = anchors_t, weight.reduction = 'cca', refdata=dR@meta.data[[cell_type_col]])
write.csv(as.data.frame(predictions_class_label), paste(output_path, 'predicted_cell_types.csv', sep='/'))

# Impute the gene expression of the query dataset
predictions_counts <- TransferData(anchorset = anchors_t, weight.reduction = 'cca', refdata=attr(attr(dR, 'assay')[[1]], 'counts'))
imputed_Q <- CreateSeuratObject(counts = predictions_counts)

# Predict the cell types in the query dataset
imputed_Q <- AddMetaData(imputed_Q, predictions_class_label['predicted.id'], col.name = 'predicted.id')
imputed_Q <- AddMetaData(imputed_Q, predictions_class_label['prediction.score.max'], col.name = 'prediction.score.max')

# Save the imputed dataset
SaveH5Seurat(imputed_Q, filename = paste(output_path, 'imputation.h5Seurat', sep='/'), overwrite = TRUE)
Convert(paste(output_path, 'imputation.h5Seurat', sep='/'), dest = "h5ad", overwrite = TRUE)
file.remove(paste(output_path, 'imputation.h5Seurat', sep='/'))

## CO-EMBED THE DATASETS AND MAKE UMAP PLOTS
# Integrate the datasets
v_ref <- c(1)
anchors <- FindIntegrationAnchors(object.list = d_list, anchor.features = features, reference=v_ref)
d_integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(d_integrated) <- "integrated"

# Co-embedding
d_integrated <- ScaleData(d_integrated)
d_integrated <- RunPCA(d_integrated)
d_integrated <- RunUMAP(d_integrated, dims = 1:50)

# Make splitted copies of the integrated data
splitted_d_integrated <- SplitObject(d_integrated, split.by='source')
integrated_query <- AddMetaData(splitted_d_integrated[["query"]], predictions_class_label['prediction.score.max'], col.name = 'prediction.score.max')
integrated_query <- AddMetaData(integrated_query, predictions_class_label['predicted.id'], col.name = 'predicted.id')
integrated_ref <- splitted_d_integrated[["reference"]]

# Plot the co-embedding
png(filename=paste(output_path, 'coembed.png', sep='/'), width=1024, height=1024)
DimPlot(d_integrated, reduction = "umap", group.by = "source")
dev.off()

# Plot the reference cell types
png(filename=paste(output_path, 'reference_cell_types.png', sep='/'), width=1024, height=1024)
DimPlot(integrated_ref, reduction = "umap", group.by = cell_type_col, label=TRUE)
dev.off()

# Plot the imputed cell types
png(filename=paste(output_path, 'imputed_cell_types.png', sep='/'), width=1024, height=1024)
DimPlot(integrated_query, reduction = "umap", group.by = 'predicted.id', label=TRUE)
dev.off()

# Plot the imputed cell type probabilities
png(filename=paste(output_path, 'imputed_cell_type_proba.png', sep='/'), width=1024, height=1024)
FeaturePlot(integrated_query, reduction='umap', features = c('prediction.score.max'), cols=c('lightgrey', 'blue'), min.cutoff='q5', max.cutoff='q95')
dev.off()

# Plot the number of detected genes in each query cell
png(filename=paste(output_path, 'query_cell_gene_counts.png', sep='/'), width=1024, height=1024)
FeaturePlot(integrated_query, reduction='umap', features = c('nFeature_RNA'), cols=c('lightgrey', 'blue'), min.cutoff='q5', max.cutoff='q95')
dev.off()

# Plot the number of detected molecules in each query cell
png(filename=paste(output_path, 'query_cell_molecule_counts.png', sep='/'), width=1024, height=1024)
FeaturePlot(integrated_query, reduction='umap', features = c('nCount_RNA'), cols=c('lightgrey', 'blue'), min.cutoff='q5', max.cutoff='q95')
dev.off()

# Plot the distributions of query data features
png(filename=paste(output_path, 'query_cell_features.png', sep='/'), width=1024, height=1024)
VlnPlot(integrated_query, c('nFeature_RNA', 'prediction.score.max'))
dev.off()