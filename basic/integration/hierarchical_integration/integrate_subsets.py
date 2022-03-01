#!/usr/bin/env python3
'''Integrate individual subsets
Local usage:
    ./integrate_subset.py integration_path cell_type_col [-d drop_gene -n n_threads]
'''

import os
from datetime import datetime
import subprocess
from optparse import OptionParser
from multiprocessing import Pool


def integrate_one_subset(subset_path, integration_script, cell_type_col, drop_gene=None, overwrite=False):
    '''Run the integration script on one subset.'''
    integrated_path = os.path.join(subset_path, 'integrated')
    done_file = os.path.join(integrated_path, 'done.txt')

    # Check if the job is already finished
    if (not overwrite) and os.path.exists(done_file):
        return

    # Run integration
    cmd = ['Rscript', '--vanilla', integration_script,
            os.path.join(subset_path, 'reference.gzip.h5seurat'),
            os.path.join(subset_path, 'query.gzip.h5seurat'),
            integrated_path, cell_type_col] 
  
    if not (drop_gene is None):
        cmd.append(drop_gene)

    print(f'Run the integration command:')
    print(' '.join(cmd))
    subprocess.check_call(cmd)

    # Mark the job as done
    with open(done_file, 'w') as f:
        f.write(f'{datetime.now()}')

def integrate_all_subsets_local(integration_path, cell_type_col, drop_gene=None, overwrite=False, n_threads=1):
    '''Run integration for all subsets on a local computer.'''
    integration_script = os.path.join(integration_path, 'integrate.R')

    # Find all subsets for integration
    integration_args = []
    for f in os.listdir(integration_path):
        subsets_path = os.path.join(integration_path, f, 'subsets')

        if os.path.exists(subsets_path):
            integration_args += [(os.path.join(subsets_path, p), 
            integration_script, cell_type_col, drop_gene, overwrite) 
            for p in os.listdir(subsets_path)]
    
    with Pool(n_threads) as p:
        p.starmap(integrate_one_subset, integration_args) 
 
def integrate_all_subsets_slurm(integration_path, cell_type_col, slurm_integration_submission_script,
        drop_gene=None, overwrite=False):
    '''Submit the integration jobs on a slurm cluster.'''
    integration_script = os.path.join(integration_path, 'integrate.R')

    if drop_gene is None:
        drop_gene = ''

    os.makedirs('slurm_job_outputs', exist_ok=True)

    # Find all subsets for integration
    integration_args = []
    for f in os.listdir(integration_path):
        subsets_path = os.path.join(integration_path, f, 'subsets')

        if not os.path.exists(subsets_path):
            continue

        for p in os.listdir(subsets_path):
            subset_path = os.path.join(subsets_path, p)
            
            # Check if the job is already finished
            done_file = os.path.join(subset_path, 'integrated', 'done.txt')
            if (not overwrite) and os.path.exists(done_file):
                continue
            
            # Submit the integration job
            if (not overwrite) and os.path.exists(done_file):
                print(f'Submit the integration job for {subset_path}.')
                cmd = ['sbatch',
                       f'--export=INTEGRATION_SCRIPT={integration_script},SUBSET_PATH={subset_path},CELL_TYPE_COL={cell_type_col},DROP_GENE={drop_gene}', 
                       slurm_integration_submission_script]
                subprocess.check_call(cmd)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-d', '--drop_gene', dest='drop_gene', action='store',
            help='The gene to drop during integration.')
    parser.add_option('-o', '--overwrite', dest='overwrite', action='store_true',
            help='Overwrite the existing result.')
    parser.add_option('-n', '--n_threads', dest='n_threads', action='store', type='int', default=1,
            help='The number of threads for a local run.')
    parser.add_option('-s', '--slurm_script', dest='slurm_script', action='store', type='string', default='',
            help='If a slurm submission script is provided, use it to run jobs on a slurm cluster.')

    (options, args) = parser.parse_args()
    integration_path, cell_type_col = args

    if len(options.slurm_script) == 0:
        integrate_all_subsets_local(integration_path, cell_type_col, drop_gene=options.drop_gene,
                overwrite=options.overwrite, n_threads=options.n_threads)

    else:
        integrate_all_subsets_slurm(integration_path, cell_type_col, options.slurm_script, drop_gene=options.drop_gene)
