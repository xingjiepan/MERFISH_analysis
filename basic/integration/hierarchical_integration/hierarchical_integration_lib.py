#!/usr/bin/env python3
'''Functions for hierarchincal_integration.
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
    os.makedirs(os.path.join(subsets_path, str(i)), exist_ok=True)
    
    output_file = os.path.join(subsets_path, str(i), f'{data_file_prefix}_{i}.gzip.h5ad')
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
    
    

