#!/usr/bin/env python3
'''Functions for prepare inputs for hierarchincal integration.
'''

import os
import subprocess
from multiprocessing import Pool

import numpy as np
import pandas as pd
import anndata
import scanpy as sc


def balanced_divide(df_cell_types, cluster_column, N_subsets, min_N_cells_per_cluster):
    '''Randomly divide the dataset into subsets.
    Arguments:
        df_cell_types: a DataFrame with the cell IDs as indices and cell types in a column.
        cluster_column: the column that stores cell types.
        N_subsets: the number of subsets to divide into.
        min_N_cells_per_cluster: the minimum number of cells for each cell type in each subset.

    Return:
        A list of list of cell IDs for each subset.
    '''
    cell_ids = [[] for i in range(N_subsets)]
    list_cluster_ids = []
    
    if cluster_column in df_cell_types.columns:
        cell_types = np.unique(df_cell_types[cluster_column])
        for ct in cell_types:
            cluster_ids = list(df_cell_types.index[df_cell_types[cluster_column] == ct])
            list_cluster_ids.append(cluster_ids)
    else:
        print(f'WARNING: The specified cluster column {cluster_column} does not exist. Continue with a single cluster.')
        list_cluster_ids.append(list(df_cell_types.index))

    for cluster_ids in list_cluster_ids:
        np.random.shuffle(cluster_ids)
        
        step_size = int(np.ceil(len(cluster_ids) / N_subsets))
        N_cells_each_selection = max(min_N_cells_per_cluster, step_size)
        
        for i in range(N_subsets):
            shift = i * step_size
            stop = min(shift + N_cells_each_selection, len(cluster_ids))
            start = max(0, stop - N_cells_each_selection)
            cell_ids[i] += cluster_ids[start : stop]
    
    return cell_ids

def generate_h5ad_to_h5seurat_conversion_R_script(output_file):
    '''Generate a R script to convert a h5ad file to a h5seurat file.'''
    script = '''# R script for converting an h5ad file to a Seurat h5 file
# Usage:
# Rscript --vanilla convert_anndata_to_seurat.R input_file

library(SeuratDisk)

args = commandArgs(trailingOnly=TRUE)
input_file = args[1]

Convert(input_file, dest='h5seurat', overwrite=TRUE, verbose=TRUE)'''

    with open(output_file, 'w') as f:
        f.write(script)

def convert_h5ad_to_h5seurat(conversion_script, input_file):
    '''Convert a h5ad file to a h5seurat file.
    NOTE: the h5ad file must be cleaned to remove non-necessary
    metadata, ortherwise the underlying R script may crash.
    '''
    conversion_cmd = ['Rscript', '--vanilla', conversion_script, input_file]
    subprocess.check_call(conversion_cmd)

def generate_one_subset(i, adata, subset_cell_ids, subsets_path, data_file_prefix, conversion_script):
    print(f'Generating the subset {data_file_prefix}_{i}.')
    adata_ds = adata[adata.obs.index.isin(subset_cell_ids[i])]
    os.makedirs(os.path.join(subsets_path, str(i), 'integrated'), exist_ok=True)
    
    output_file = os.path.join(subsets_path, str(i), f'{data_file_prefix}.gzip.h5ad')
    adata_ds.write(output_file, compression='gzip')
    convert_h5ad_to_h5seurat(conversion_script, output_file)
    os.remove(output_file)

def generate_subsets(adata, N_subsets, n_repeat, N_subsets_to_write, 
        cell_type_column, min_N_cells_per_cluster,
        subsets_path, conversion_script, data_file_prefix,
        n_threads=1):
    '''Generate the subsets of a dataset.'''
    
    print(f'Determining the cell IDs for {data_file_prefix} subsets.')
    subset_cell_ids = []
    for i in range(n_repeat):
        subset_cell_ids += balanced_divide(adata.obs, cell_type_column, N_subsets, min_N_cells_per_cluster)

    # Generate the subsets in prallel
    with Pool(n_threads) as p:
        p.starmap(generate_one_subset, [(i, adata, subset_cell_ids, subsets_path, data_file_prefix, conversion_script) 
            for i in range(N_subsets_to_write)])

def prepare_integration_inputs(output_path, reference_adata_file, query_adata_file, 
        reference_cell_type_column, query_cell_type_column, approximate_subset_size=10000,
        n_repeat_query=3, min_N_cells_per_cluster=50, n_threads=1):
    '''Prepare the inputs for the actual integration script.
    NOTE: the adata files should already be cleaned to remove all
    unnecessary metadata.
    '''
    # Load the data
    reference_adata = sc.read_h5ad(reference_adata_file)
    query_adata = sc.read_h5ad(query_adata_file)

    # Determine the number of subsets to generate
    N_subsets_reference = max(1, int(reference_adata.shape[0] / approximate_subset_size))
    N_subsets_query = max(1, int(query_adata.shape[0] / approximate_subset_size))
    
    n_repeat_reference = int(np.ceil(N_subsets_query * n_repeat_query / N_subsets_reference))
    N_subsets_to_write = N_subsets_query * n_repeat_query
    print(f'Generating {N_subsets_to_write} subsets in total.')

    # Generate the subsets
    subsets_path = os.path.join(output_path, 'subsets')
    os.makedirs(subsets_path, exist_ok=True)
    conversion_script = os.path.join(output_path, 'convert_anndata_to_seurat.R') 
    generate_h5ad_to_h5seurat_conversion_R_script(conversion_script)

    generate_subsets(reference_adata, N_subsets_reference, n_repeat_reference, N_subsets_to_write,
            reference_cell_type_column, min_N_cells_per_cluster, subsets_path, conversion_script, 'reference',
            n_threads=n_threads)

    generate_subsets(query_adata, N_subsets_query, n_repeat_query, N_subsets_to_write,
            query_cell_type_column, min_N_cells_per_cluster, subsets_path, conversion_script, 'query',
            n_threads=n_threads)
    
def generate_seurat_integration_script(output_file, impute_gene_expression=True, plot_coembedding=True):
    '''Generate an R script for integration using Seurat.'''
    script = '''# R script for integrating a query Seurat dataset to a reference Seurat dataset.
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

# Normalize and identify variable features for each dataset independently
d_list <- lapply(X = d_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = d_list)
features = features[features != drop_gene]

## IMPUTATION AND LABEL TRANSFER
# Find the anchors for imputation and label transfer
anchors_t <- FindTransferAnchors(reference = dR, query = dQ, reduction='cca', features=features)

# Predict the class labels and same the predictions as a csv file 
predictions_class_label <- TransferData(anchorset = anchors_t, weight.reduction = 'cca', refdata=dR@meta.data[[cell_type_col]])
write.csv(as.data.frame(predictions_class_label), paste(output_path, 'predicted_cell_types.csv', sep='/'))
'''
    if impute_gene_expression:
        script += '''
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
'''
    if plot_coembedding:
        script += '''
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
dev.off()'''

    with open(output_file, 'w') as f:
        f.write(script)
