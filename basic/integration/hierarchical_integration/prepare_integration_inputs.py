#!/usr/bin/env python3
'''Generate the input files for hierarchical integration.
'''
import os
import subprocess
from multiprocessing import Pool
from optparse import OptionParser

import numpy as np
import pandas as pd
import scipy.sparse
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

def generate_h5seurat_conversion_R_script(script_file):
    '''Generate a R script to convert a h5ad file to a h5seurat file.'''
    script = '''# R script for converting an h5ad file to a Seurat h5 file
# Usage:
# Rscript --vanilla convert_anndata_to_seurat.R mtx_file cells_file features_file metadata_file output_file

library(Seurat)
library(SeuratDisk)
library(Matrix)
library(reticulate)

args = commandArgs(trailingOnly=TRUE)
mtx_file = args[1]
cells_file = args[2]
features_file = args[3]
metadata_file = args[4]
output_file = args[5]

cells = read.csv(cells_file, colClasses=c("character"), header=FALSE)
features = read.csv(features_file, colClasses=c("character"), header=FALSE)
metadata = read.csv(metadata_file)
attr(metadata, "row.names") = cells[[1]]

scipy = import("scipy.sparse")
m = scipy$load_npz(mtx_file)

m@Dimnames[[1]] = cells[[1]]
m@Dimnames[[2]] = features[[1]]

seurat_object <- CreateSeuratObject(counts=t(m), meta.data=metadata)

SaveH5Seurat(seurat_object, output_file, overwrite=TRUE)'''

    with open(script_file, 'w') as f:
        f.write(script)

def convert_anndata_to_h5seurat_file(adata, conversion_script, scratch_path, output_file):
    '''Convert and anndata to a h5seurat file.'''
    # Generate the intermediate files
    mtx_file = os.path.join(scratch_path, 'counts.npz') 
    cells_file = os.path.join(scratch_path, 'cells.csv') 
    features_file = os.path.join(scratch_path, 'features.csv') 
    metadata_file = os.path.join(scratch_path, 'metadata.csv') 

    print('Write the matrix file.')
    scipy.sparse.save_npz(mtx_file, scipy.sparse.csc_matrix(adata.X))
    print('Write the cells file.')
    adata.obs[[]].to_csv(cells_file, header=False)
    print('Write the features file.')
    adata.var[[]].to_csv(features_file, header=False)
    print('Write the metadata file.')
    adata.obs.to_csv(metadata_file, index=False)

    # Run the conversion script
    conversion_cmd = ['Rscript', '--vanilla', conversion_script, 
            mtx_file, cells_file, features_file, metadata_file, output_file]
    print('Running command:')
    print(' '.join(conversion_cmd))
    subprocess.check_call(conversion_cmd)

    # Remove the intermediate files
    os.remove(mtx_file)
    os.remove(cells_file)
    os.remove(features_file)
    os.remove(metadata_file)

def all_cell_id_files_already_exist(N_subsets_to_write, subsets_path, data_file_prefix):
    all_exist = True

    for i in range(N_subsets_to_write):
        cell_id_file = os.path.join(subsets_path, str(i), f'{data_file_prefix}_ids.csv')
        if not os.path.exists(cell_id_file): 
            all_exist = False

    if all_exist:
        print(f'All cell id files under {subsets_path} already exist.')
    return all_exist

def save_cells_ids_for_subsets(i, subset_cell_ids, subsets_path, data_file_prefix, overwrite=True):
    '''Save the cell ids for each subset of integration.'''
    os.makedirs(os.path.join(subsets_path, str(i), 'integrated'), exist_ok=True)
    cell_id_file = os.path.join(subsets_path, str(i), f'{data_file_prefix}_ids.csv')

    if os.path.exists(cell_id_file) and (not overwrite):
        return

    pd.DataFrame({'cell_ids':subset_cell_ids}).to_csv(cell_id_file, header=False, index=False)

def generate_one_subset(i, adata, subsets_path, data_file_prefix, conversion_script, overwrite=False):
    final_output = os.path.join(subsets_path, str(i), f'{data_file_prefix}.gzip.h5seurat')
    if os.path.exists(final_output) and (not overwrite):
        print(f'The file {final_output} already exists.')
        return

    print(f'Generating the subset {data_file_prefix}_{i}.')
    cell_id_file = os.path.join(subsets_path, str(i), f'{data_file_prefix}_ids.csv')
    subset_cell_ids =  np.array(pd.read_csv(cell_id_file, header=None).iloc[:, 0])

    adata_ds = adata[adata.obs.index.isin(subset_cell_ids)]
    os.makedirs(os.path.join(subsets_path, str(i), 'integrated'), exist_ok=True)
    
    scratch_path = os.path.join(subsets_path, str(i))
    output_file = os.path.join(subsets_path, str(i), f'{data_file_prefix}.gzip.h5seurat')
    convert_anndata_to_h5seurat_file(adata_ds, conversion_script, scratch_path, output_file)

def generate_subsets(adata, N_subsets, n_repeat, N_subsets_to_write, 
        cell_type_column, min_N_cells_per_cluster,
        subsets_path, conversion_script, data_file_prefix,
        n_threads=1, overwrite=False):
    '''Generate the subsets of a dataset.'''
    
    print(f'Determining the cell IDs for {data_file_prefix} subsets.')
    # Check if cell IDs for each subset are already generated
    all_exist = all_cell_id_files_already_exist(N_subsets_to_write, subsets_path, data_file_prefix) 

    # If not all cell_id files already exist, generate all files 
    if not all_exist:
        all_subset_cell_ids = []
        for i in range(n_repeat):
            all_subset_cell_ids += balanced_divide(adata.obs, cell_type_column, N_subsets, min_N_cells_per_cluster)

        # Save the cell IDs for each subset
        for i in range(N_subsets_to_write):
            save_cells_ids_for_subsets(i, all_subset_cell_ids[i], subsets_path, data_file_prefix)

    # Generate the subsets sequentially
    if n_threads == 1:
        for i in range(N_subsets_to_write):
            generate_one_subset(i, adata, subsets_path, data_file_prefix, conversion_script)

    # Generate the subsets in prallel
    else:
        with Pool(n_threads) as p:
            p.starmap(generate_one_subset, [(i, adata, subsets_path, data_file_prefix, conversion_script) 
                for i in range(N_subsets_to_write)])

def prepare_integration_inputs_for_one_cell_type(output_path, reference_adata, query_adata, 
        reference_cell_type_column, query_cell_type_column, approximate_subset_size=10000,
        n_repeat_query=3, min_N_cells_per_cluster=50, n_threads=1):
    '''Prepare the inputs for the actual integration script.
    NOTE: the adata files should already be cleaned to remove all
    unnecessary metadata.
    '''
    # Adjust the approximate_subset_size such that the query and reference datasets have similar sizes
    approximate_subset_size = min([approximate_subset_size, reference_adata.shape[0], query_adata.shape[0]])

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
    generate_h5seurat_conversion_R_script(conversion_script)

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
    for i, ct in enumerate(cell_types_to_split): 
        print(f'\nGenerating integration inputs for {reference_col_to_split}:{ct}.')
        print(f'This is the {i+1} / {len(cell_types_to_split)} cell types.') 

        output_path_ct = os.path.join(output_path, ct)
        reference_adata_ct = reference_adata[reference_adata.obs[reference_col_to_split] == ct]
        query_adata_ct = query_adata[query_adata.obs[query_col_to_split] == ct]

        # Only proceed if the cell type exists in the query and reference
        if reference_adata_ct.shape[0] == 0 or query_adata_ct.shape[0] == 0:
            continue

        prepare_integration_inputs_for_one_cell_type(output_path_ct, reference_adata_ct, query_adata_ct,
                reference_cell_type_column, query_cell_type_column, approximate_subset_size=approximate_subset_size,
                n_repeat_query=n_repeat_query, min_N_cells_per_cluster=min_N_cells_per_cluster, n_threads=n_threads)

def generate_seurat_integration_script(output_file, impute_gene_expression=True, plot_coembedding=True,
                                       continuous_columns_to_impute=[], variable_genes_fraction=0.5):
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


## DEFINE THE GENES FOR INTEGRATION
# Only keep the genes in the query dataset for normalization
all_query_genes <- rownames(dQ)
de_gene_file = paste(output_path, 'preselected_de_genes.csv', sep='/') 

if (file.exists(de_gene_file)){ # Use preselected DE genes
    print('Load preselected differentially expressed genes.')
    de_df = read.csv(de_gene_file)
    features = intersect(de_df$de_genes, all_query_genes)
     
}else{ # Find DE genes.
    print('Look for differentially expressed genes.')
    dR_query_genes <- subset(dR, features = all_query_genes)

    d_list <- list(dR_query_genes, dQ)

    # Normalize and identify variable features for each dataset independently
    d_list <- lapply(X = d_list, FUN = function(x) {
        x <- NormalizeData(x)
'''
    script += f'''        x <- FindVariableFeatures(x, selection.method="vst", nfeatures=as.integer({variable_genes_fraction} * length(all_query_genes)))
'''

    script += '''   })

    # Select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = d_list)
    features = features[features != drop_gene]
}

print('Integration using the following genes:')
print(features)

# Create Seurat objects for integration 
dR_de_genes = subset(dR, features = features)
dQ_de_genes = subset(dQ, features = features)
d_list <- list(dR_de_genes, dQ_de_genes)
d_list <- lapply(X = d_list, FUN = function(x) {
  x <- NormalizeData(x)
})

## MAKE THE CO-EMBEDDING PLOT BEFORE INTEGRATION
d_merged <- merge(dR_de_genes, y = dQ_de_genes, add.cell.ids = c("reference", "query"))
d_merged <- ScaleData(d_merged)
d_merged <- RunPCA(d_merged, features=features)
d_merged <- RunUMAP(d_merged, dims = 1:min(length(d_merged@reductions[["pca"]]), 50))

png(filename=paste(output_path, 'coembed_before_integration.png', sep='/'), width=1024, height=1024)
DimPlot(d_merged, reduction = "umap", group.by = "source")
dev.off()

## IMPUTATION AND LABEL TRANSFER
# Find the anchors for imputation and label transfer
anchors_t <- FindTransferAnchors(reference = dR_de_genes, query = dQ_de_genes, reduction='cca', features=features, dims=1:min(length(features) - 1, 30))

# Predict the class labels and same the predictions as a csv file 
predictions_class_label <- TransferData(anchorset = anchors_t, weight.reduction = 'cca', refdata=as.factor(dR@meta.data[[cell_type_col]]))
write.csv(as.data.frame(predictions_class_label), paste(output_path, 'predicted_cell_types.csv', sep='/'))
'''
    if len(continuous_columns_to_impute) > 0:
        cont_col_string = 'c(' + ','.join(["'" + c + "'" for c in continuous_columns_to_impute]) + ')'
        script += f'''
# Impute the continuous values in selected columns
ref_cont_values = t(data.matrix(dR@meta.data[{cont_col_string}]))
predicted_cont_values = TransferData(anchorset=anchors_t, weight.reduction='cca', refdata=ref_cont_values)
write.csv(as.data.frame(Transpose(predicted_cont_values@data)), paste(output_path, 'predicted_cont_values.csv', sep='/'))
'''
    if impute_gene_expression:
        script += '''
# Impute the gene expression of the query dataset
predictions_counts <- TransferData(anchorset = anchors_t, weight.reduction = 'cca', refdata=attr(attr(dR, 'assay')[[1]], 'counts'))
imputed_Q <- CreateSeuratObject(counts = predictions_counts)

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
anchors <- FindIntegrationAnchors(object.list = d_list, anchor.features = features, reference=v_ref,
                                 dims=1:min(length(features) - 1, 30))
d_integrated <- IntegrateData(anchorset = anchors, dims=1:min(length(features) - 1, 30))
DefaultAssay(d_integrated) <- "integrated"

# Co-embedding
d_integrated <- ScaleData(d_integrated)
d_integrated <- RunPCA(d_integrated)
d_integrated <- RunUMAP(d_integrated, dims = 1:min(length(d_integrated@reductions[["pca"]]), 50))

# Calculate and save the mixing score for the integration
d_integrated <- AddMetaData(d_integrated, MixingMetric(d_integrated, 'source', reduction='pca', 
                            dims=1:min(length(d_integrated@reductions[["pca"]]), 20)), 'mixing_score')

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
    parser.add_option('--continuous_columns_to_impute', dest='continuous_columns_to_impute', 
            action='store', type='string', default='',
            help='A comma separated list of column names for continuous variables to be imputed')

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
    if len(options.continuous_columns_to_impute) > 0:
        continuous_columns_to_impute = options.continuous_columns_to_impute.split(',')
    else: 
        continuous_columns_to_impute = []
    
    generate_seurat_integration_script(os.path.join(output_path, 'integrate.R'), 
            impute_gene_expression=impute_gene_expression, plot_coembedding=True,
            continuous_columns_to_impute=continuous_columns_to_impute)

    # Generate the input data
    prepare_integration_inputs_for_one_round(output_path, reference_adata_file, query_adata_file, 
        reference_col_to_split, query_col_to_split, 
        reference_col_cell_type, None, approximate_subset_size=approximate_subset_size,
        n_repeat_query=n_repeat_query, min_N_cells_per_cluster=min_N_cells_per_cluster, n_threads=n_threads)

    print('Finished preparing inputs.')
    with open(os.path.join(output_path, 'prepare_inputs.done'), 'w') as f:
        f.write('')
