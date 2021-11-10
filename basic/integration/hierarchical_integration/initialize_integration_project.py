#!/usr/bin/env python3

import os

import numpy as np
import anndata



def load_and_clean_adata(adata_file, obs_columns_to_keep):
    '''Clean an adata file by removing all metadata except for specified
    columns.'''
    adata = anndata.read_h5ad(adata_file)
    adata_cleaned = anndata.AnnData(X=adata.X, obs=adata.obs[obs_columns_to_keep], var=adata.var[[]])
    return adata_cleaned

def generate_script_for_prepare_inputs(script_home, project_path, rd, 
        reference_adata_file, query_adata_file, 
        reference_col_to_split, query_col_to_split, reference_col_cell_type,
        approximate_subset_size, n_repeat_query=3, min_N_cells_per_cluster=30, n_threads=1,
        impute_gene_expression=False):
    '''Generate script for preparing inputs of a round.'''
    
    task_script = os.path.abspath(os.path.join(script_home, 'prepare_integration_inputs.py'))
    output_path = os.path.abspath(os.path.join(project_path, f'round{rd}'))
    reference_adata_file = os.path.abspath(reference_adata_file)
    query_adata_file = os.path.abspath(query_adata_file)

    cmd = [task_script, '-r', str(n_repeat_query), '-m', str(min_N_cells_per_cluster), 
            '-n', str(n_threads)]
    if impute_gene_expression:
        cmd.appen('-i')

    cmd += [output_path, reference_adata_file, query_adata_file, reference_col_to_split,
            query_col_to_split, reference_col_cell_type, str(approximate_subset_size)]

    script = f'''#!/bin/bash
    
{' '.join(cmd)}
    '''

    with open(os.path.join(project_path, f'round{rd}_prepare_inputs.sh'), 'w') as f:
        f.write(script)

def generate_script_for_integrate_subsets(script_home, project_path, rd, reference_col_cell_type, 
        drop_gene=None, overwrite=False, n_threads=1, slurm=False):
    '''Generate script for integrating subsets of a round.'''
    task_script = os.path.abspath(os.path.join(script_home, 'integrate_subsets.py'))
    round_path = os.path.abspath(os.path.join(project_path, f'round{rd}'))

    cmd = [task_script, '-n', str(n_threads)]
    
    if not (drop_gene is None):
        cmd += ['-d', drop_gene]
    if overwrite:
        cmd.append('-o')
    if slurm:
        cmd.append('-s')

    cmd += [round_path, reference_col_cell_type]

    script = f'''#!/bin/bash
    
{' '.join(cmd)}
    '''

    with open(os.path.join(project_path, f'round{rd}_integrate_subsets.sh'), 'w') as f:
        f.write(script)

def generate_script_for_analyze_result(script_home, project_path, rd, query_adata_file, new_query_adata_file,
        prediction_col, prediction_proba_col):
    '''Generate script for analyzing result of a round.'''
    pass

def initialize_integration_project(script_home, project_path, reference_adata_file, query_adata_file,
        reference_columns_by_rounds):
    '''Initialize an integration project.'''
    # Create the working directory and copy the cleaned anndata files here
    os.makedirs(project_path, exist_ok=True)

    # Save the cleaned data to the project directory
    cleaned_reference_adata_file = os.path.join(project_path, 'reference.h5ad')
    reference_adata_cleaned = load_and_clean_adata(reference_adata_file, reference_columns_by_rounds)
    reference_adata_cleaned.obs['root_type'] = 'all'
    reference_adata_cleaned.write(cleaned_reference_adata_file)

    cleaned_query_adata_file = os.path.join(project_path, 'query.h5ad')
    query_adata_cleaned = load_and_clean_adata(query_adata_file, [])
    query_adata_cleaned.obs['root_type'] = 'all'
    query_adata_cleaned.write(cleaned_query_adata_file)

    # Create the directory structure for each round of integration

    for i, col in enumerate(reference_columns_by_rounds):
        round_dir = os.path.join(project_path, f'round{i}')
        os.makedirs(round_dir, exist_ok=True)

    # Generate scripts for integration
    
    for i, col in enumerate(reference_columns_by_rounds):

        # Generate input preparation scripts
        if i == 0:
            query_adata_file = cleaned_query_adata_file
            reference_col_to_split = 'root_type'
            query_col_to_split = 'root_type'
       
        else:
            query_adata_file = os.path.join(project_path, f'round{i-1}', 'integrated.h5ad')
            reference_col_to_split = reference_columns_by_rounds[i - 1]
            query_col_to_split = 'prediction_' + reference_col_to_split

        reference_adata_file = cleaned_reference_adata_file
        reference_col_cell_type = col
        approximate_subset_size = 1000 #TODO
        n_repeat_query=3 #TODO
        min_N_cells_per_cluster=30 #TODO
        n_threads=16 #TODO
        impute_gene_expression=False #TODO

        generate_script_for_prepare_inputs(script_home, project_path, i, 
            reference_adata_file, query_adata_file, 
            reference_col_to_split, query_col_to_split, reference_col_cell_type,
            approximate_subset_size, n_repeat_query, min_N_cells_per_cluster, n_threads,
            impute_gene_expression)

        # Generate actual integration script
        
        drop_gene=None #TODO
        overwrite=False #TODO
        n_threads=16 #TODO
        slurm=False #TODO
        generate_script_for_integrate_subsets(script_home, project_path, i, reference_col_cell_type, 
                drop_gene=drop_gene, overwrite=overwrite, n_threads=n_threads, slurm=slurm)
        


if __name__ == '__main__':
    script_home = '.'
    project_path = 'test/test_integration_project'
    reference_adata_file = 'test/scRNAseq_downsample_0.gzip.h5ad'
    query_adata_file = 'test/merfish_downsample_0.gzip.h5ad'

    reference_columns_by_rounds = ['cell_class1', 'seurat_clusters']
    

    initialize_integration_project(script_home, project_path, reference_adata_file, 
            query_adata_file, reference_columns_by_rounds)

