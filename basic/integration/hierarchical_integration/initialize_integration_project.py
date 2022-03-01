#!/usr/bin/env python3

import os
from optparse import OptionParser

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
        impute_gene_expression=False, continuous_columns_to_impute=[], harvard_rc=False):
    '''Generate script for preparing inputs of a round.'''
    
    task_script = os.path.abspath(os.path.join(script_home, 'prepare_integration_inputs.py'))
    output_path = os.path.abspath(os.path.join(project_path, f'round{rd}'))
    reference_adata_file = os.path.abspath(reference_adata_file)
    query_adata_file = os.path.abspath(query_adata_file)

    cmd = [task_script, '-r', str(n_repeat_query), '-m', str(min_N_cells_per_cluster), 
            '-n', str(n_threads)]
    if impute_gene_expression:
        cmd.append('-i')
    
    if len(continuous_columns_to_impute) > 0:
        cmd += ['--continuous_columns_to_impute', ','.join(continuous_columns_to_impute)]

    cmd += [output_path, reference_adata_file, query_adata_file, reference_col_to_split,
            query_col_to_split, reference_col_cell_type, str(approximate_subset_size)]

    if harvard_rc:
        script = f'''#!/bin/bash
#SBATCH --job-name=integrate_
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-24:00          # Runtime in D-HH:MM
#SBATCH -p zhuang   # Partition to submit to
#SBATCH --mem-per-cpu=48000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o slurm_job_outputs/output_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e slurm_job_outputs/output_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --account=zhuang_lab
#SBATCH --exclude=holyzhuang01,holy2a15313,holy2c[093301,093302,093401,093402,01213,01214] # Some nodes fail to run the job

# Load the module for the correct version of hdf5
module load gcc/10.2.0-fasrc01 openmpi/4.1.1-fasrc01 hdf5/1.12.1-fasrc01

# Load the module for R
module load R_core/4.0.5-fasrc01

# Load the R packages
module load R_packages/4.0.5-fasrc01

{' '.join(cmd)}
'''
    os.makedirs(os.path.join(project_path, 'slurm_job_outputs')), exist_ok=True)
    else:
        script = f'''#!/bin/bash
    
{' '.join(cmd)}
    '''

    with open(os.path.join(project_path, f'round{rd}_prepare_inputs.sh'), 'w') as f:
        f.write(script)

def generate_script_for_integrate_subsets(script_home, project_path, rd, reference_col_cell_type, 
        drop_gene=None, overwrite=False, n_threads=1, slurm_script=''):
    '''Generate script for integrating subsets of a round.'''
    task_script = os.path.abspath(os.path.join(script_home, 'integrate_subsets.py'))
    round_path = os.path.abspath(os.path.join(project_path, f'round{rd}'))

    cmd = [task_script, '-n', str(n_threads)]
    
    if not (drop_gene is None):
        cmd += ['-d', drop_gene]
    if overwrite:
        cmd.append('-o')
    if len(slurm_script) > 0:
        cmd += ['-s', os.path.abspath(slurm_script)]

    cmd += [round_path, reference_col_cell_type]

    script = f'''#!/bin/bash
    
{' '.join(cmd)}
    '''

    with open(os.path.join(project_path, f'round{rd}_integrate_subsets.sh'), 'w') as f:
        f.write(script)

def generate_script_for_analyze_result(script_home, project_path, rd,
        adata_file_before_integration, prediction_cell_type_col, confidence_threshold):
    '''Generate script for analyzing result of a round.'''
    task_script = os.path.abspath(os.path.join(script_home, 'analyze_integratioin_results.py'))
    round_path = os.path.abspath(os.path.join(project_path, f'round{rd}'))
    adata_file_before_integration = os.path.abspath(adata_file_before_integration)
   
    cmd = [task_script, '-t', str(confidence_threshold),
            round_path, adata_file_before_integration, prediction_cell_type_col]

    script = f'''#!/bin/bash
    
{' '.join(cmd)}
    '''

    with open(os.path.join(project_path, f'round{rd}_analyze_integration_result.sh'), 'w') as f:
        f.write(script)

def initialize_integration_project(options):
    '''Initialize an integration project.'''
    # Parse the reference columns
    reference_columns_by_rounds = options.reference_columns_by_rounds.split(',')
    
    if not (options.continuous_columns_to_impute is None):
        continuous_columns_to_impute = options.continuous_columns_to_impute.split(',')
    else: 
        continuous_columns_to_impute = []

    assert(not ('root_type' in reference_columns_by_rounds))

    # Create the working directory and copy the cleaned anndata files here
    os.makedirs(options.project_path, exist_ok=True)

    # Save the cleaned data to the project directory
    cleaned_reference_adata_file = os.path.join(options.project_path, 'reference.h5ad')
    reference_adata_cleaned = load_and_clean_adata(options.reference_adata_file, 
            reference_columns_by_rounds + continuous_columns_to_impute)
    reference_adata_cleaned.obs['root_type'] = 'all'
    reference_adata_cleaned.write(cleaned_reference_adata_file)

    cleaned_query_adata_file = os.path.join(options.project_path, 'query.h5ad')
    query_adata_cleaned = load_and_clean_adata(options.query_adata_file, [])
    query_adata_cleaned.obs['root_type'] = 'all'
    query_adata_cleaned.write(cleaned_query_adata_file)

    # Create the directory structure for each round of integration
    for i, col in enumerate(reference_columns_by_rounds):
        round_dir = os.path.join(options.project_path, f'round{i}')
        os.makedirs(round_dir, exist_ok=True)

    # Generate scripts for integration
    for i, col in enumerate(reference_columns_by_rounds):

        impute_gene_expression = False # TODO: implement gene expression imputation

        # The first round has some different parameters
        if i == 0:
            query_adata_file_current = cleaned_query_adata_file
            reference_col_to_split = 'root_type'
            query_col_to_split = 'root_type'
       
        else:
            query_adata_file_current = os.path.join(options.project_path, f'round{i-1}', 'integrated.h5ad')
            reference_col_to_split = reference_columns_by_rounds[i - 1]
            query_col_to_split = 'prediction_' + reference_col_to_split + '_filtered'

        # Generate input preparation scripts
        if len(options.slurm_script) > 0:
            harvard_rc = True
        else:
            harvard_rc = False

        generate_script_for_prepare_inputs(options.script_home, options.project_path, i, 
            cleaned_reference_adata_file, query_adata_file_current, 
            reference_col_to_split, query_col_to_split, col, options.approximate_subset_size, 
            options.n_repeat_query, options.min_N_cells_per_cluster, options.n_threads,
            impute_gene_expression, continuous_columns_to_impute=continuous_columns_to_impute,
            harvard_rc=harvard_rc)

        # Generate actual integration script
        generate_script_for_integrate_subsets(options.script_home, options.project_path, i, 
                col, drop_gene=options.drop_gene, overwrite=options.overwrite, 
                n_threads=options.n_threads, slurm_script=options.slurm_script)

        # Generate the integration analysis script
        generate_script_for_analyze_result(options.script_home, options.project_path, i,
            query_adata_file_current, 'prediction_' + col, options.mixing_threshold)
        


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-i', '--script_home', dest='script_home', action='store', type='string',
            help='The directory of the integration scripts.')
    parser.add_option('-p', '--project_path', dest='project_path', action='store', type='string',
            help='The path to the integration project.')
    parser.add_option('-q', '--query_adata_file', dest='query_adata_file', action='store', type='string',
            help='The query h5ad file.')
    parser.add_option('-r', '--reference_adata_file', dest='reference_adata_file', action='store', type='string',
            help='The reference h5ad file.')
    parser.add_option('-c', '--reference_columns_by_rounds', dest='reference_columns_by_rounds', 
            action='store', type='string', default='',
            help='A comma separated list of column names for each round of integration.')
    parser.add_option('--continuous_columns_to_impute', dest='continuous_columns_to_impute', action='store', type='string',
            help='A comma separated list of column names for continuous variables to be imputed')
    parser.add_option('-a', '--approximate_subset_size', dest='approximate_subset_size', action='store', type='int', default=10000,
            help='The approximate size of each subset for integration.')
    parser.add_option('-e', '--n_repeat_query', dest='n_repeat_query', action='store', type='int', default=3,
            help='The number of repeat for each cell in the query dataset.')
    parser.add_option('-m', '--min_N_cells_per_cluster', dest='min_N_cells_per_cluster', action='store', type='int', default=50,
            help='The minimum number of cells in each downsampled subsets.')
    parser.add_option('-t', '--mixing_threshold', dest='mixing_threshold', action='store', type='float',
            help='The mixing score threshold for filtering integrated cells.')
    parser.add_option('-d', '--drop_gene', dest='drop_gene', action='store',
            help='The gene to drop during integration.')
    parser.add_option('-o', '--overwrite', dest='overwrite', action='store_true',
            help='Overwrite the existing result.')
    parser.add_option('-s', '--slurm_script', dest='slurm_script', action='store', type='string', default='',
            help='If a slurm submission script is provided, use it to run jobs on a slurm cluster.')
    parser.add_option('-n', '--n_threads', dest='n_threads', action='store', type='int', default=1,
            help='The number of threads for a local run.')

    (options, args) = parser.parse_args()

    initialize_integration_project(options)

