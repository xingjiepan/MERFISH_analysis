#!/bin/bash
#SBATCH --job-name=integrate_
#SBATCH -c 1                # Number of cores (-c)
#SBATCH -t 0-03:00          # Runtime in D-HH:MM
#SBATCH -p zhuang   # Partition to submit to
#SBATCH --mem-per-cpu=10000           # Memory pool for all cores (see also --mem-per-cpu)
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

Rscript --vanilla ${INTEGRATION_SCRIPT} ${SUBSET_PATH}/reference.gzip.h5seurat ${SUBSET_PATH}/query.gzip.h5seurat ${SUBSET_PATH}/integrated ${CELL_TYPE_COL} ${DROP_GENE}
