#!/usr/bin/env python3
'''Generate the input files for hierarchical integration.
'''
import os
import subprocess
from multiprocessing import Pool
from optparse import OptionParser

import numpy as np
import pandas as pd
import anndata


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
    print('Running command:')
    print(' '.join(conversion_cmd))
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

def prepare_integration_inputs_for_one_cell_type(output_path, reference_adata, query_adata, 
        reference_cell_type_column, query_cell_type_column, approximate_subset_size=10000,
        n_repeat_query=3, min_N_cells_per_cluster=50, n_threads=1):
    '''Prepare the inputs for the actual integration script.
    NOTE: the adata files should already be cleaned to remove all
    unnecessary metadata.
    '''
    # Determine the number of subsets to generate
    N_subsets_reference = max(1, int(np.ceil(reference_adata.shape[0] / approximate_subset_size)))
    N_subsets_query = max(1, int(np.ceil(query_adata.shape[0] / approximate_subset_size)))
    
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
   
def prepare_integration_inputs_for_one_round(output_path, reference_adata_file, query_adata_file, 
        reference_col_to_split, query_col_to_split, 
        reference_cell_type_column, query_cell_type_column, approximate_subset_size=10000,
        n_repeat_query=3, min_N_cells_per_cluster=50, n_threads=1):
    '''Prepare inputs for one round of integration.'''
    # Load the data
    reference_adata = anndata.read_h5ad(reference_adata_file)
    query_adata = anndata.read_h5ad(query_adata_file)

    # Get the cell types by which the dataset is splitted into smaller sets
    cell_types_to_split = np.unique(reference_adata.obs[reference_col_to_split])

    # Generate inputs for each splitted cell type
    for ct in cell_types_to_split: 
        print(f'Generating integration inputs for {reference_col_to_split}:{ct}.')
        
        output_path_ct = os.path.join(output_path, ct)
        reference_adata_ct = reference_adata[reference_adata.obs[reference_col_to_split] == ct]
        query_adata_ct = query_adata[query_adata.obs[query_col_to_split] == ct]
        #TODO: deal with the case when there is no cell in the query

        prepare_integration_inputs_for_one_cell_type(output_path_ct, reference_adata_ct, query_adata_ct,
                reference_cell_type_column, query_cell_type_column, approximate_subset_size=approximate_subset_size,
                n_repeat_query=n_repeat_query, min_N_cells_per_cluster=min_N_cells_per_cluster, n_threads=n_threads)

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

# Calculate and save the mixing score for the integration
d_integrated <- AddMetaData(d_integrated, MixingMetric(d_integrated, 'source', reduction='pca', dims=1:20), 'mixing_score')

# Make splitted copies of the integrated data
splitted_d_integrated <- SplitObject(d_integrated, split.by='source')
integrated_query <- AddMetaData(splitted_d_integrated[["query"]], predictions_class_label['prediction.score.max'], col.name = 'prediction.score.max')
integrated_query <- AddMetaData(integrated_query, predictions_class_label['predicted.id'], col.name = 'predicted.id')
integrated_ref <- splitted_d_integrated[["reference"]]

# Save the mixing scores of the query cells
write.csv(integrated_query@meta.data["mixing_score"], paste(output_path, 'query_cell_mixing_scores.csv', sep='/'))

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

# Plot the mixing score each query cell
png(filename=paste(output_path, 'query_cell_mixing_scores.png', sep='/'), width=1024, height=1024)
FeaturePlot(integrated_query, reduction='umap', features = c('mixing_score'), cols=c('lightgrey', 'blue'), min.cutoff='q5', max.cutoff='q95')
dev.off()

# Plot the distributions of query data features
png(filename=paste(output_path, 'query_cell_features.png', sep='/'), width=1024, height=1024)
VlnPlot(integrated_query, c('nFeature_RNA', 'prediction.score.max', 'mixing_score'))
dev.off()'''

    with open(output_file, 'w') as f:
        f.write(script)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-r', '--n_repeat_query', dest='n_repeat_query', action='store', type='int', default=3,
            help='Number of repeats for each query cell.')
    parser.add_option('-m', '--min_N_cells_per_cluster', dest='min_N_cells_per_cluster', action='store', 
            type='int', default=30, help='Minimal number of cells that should exist in each reference cluster.')
    parser.add_option('-n', '--n_threads', dest='n_threads', action='store', type='int', default=1,
            help='The gene to drop during integration.')
    parser.add_option('-i', '--impute', dest='impute', action='store_true',
            help='Impute gene expression.')

    (options, args) = parser.parse_args()
    
    output_path = args[0] 
    reference_adata_file = args[1] 
    query_adata_file = args[2]
    reference_col_to_split = args[3]
    query_col_to_split = args[4]
    reference_col_cell_type = args[5]
    approximate_subset_size = int(args[6])

    n_repeat_query = options.n_repeat_query
    min_N_cells_per_cluster = options.min_N_cells_per_cluster
    n_threads = options.n_threads
    impute_gene_expression = options.impute

    # Generate the integration script
    generate_seurat_integration_script(os.path.join(output_path, 'integrate.R'), 
            impute_gene_expression=impute_gene_expression, plot_coembedding=True)

    # Generate the input data
    prepare_integration_inputs_for_one_round(output_path, reference_adata_file, query_adata_file, 
        reference_col_to_split, query_col_to_split, 
        reference_col_cell_type, None, approximate_subset_size=approximate_subset_size,
        n_repeat_query=n_repeat_query, min_N_cells_per_cluster=min_N_cells_per_cluster, n_threads=n_threads)

