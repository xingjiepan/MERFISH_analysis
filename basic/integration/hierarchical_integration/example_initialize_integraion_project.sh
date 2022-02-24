#!/bin/bash

./initialize_integration_project.py -i . -p test/test_integration_project -q test/merfish_downsample_0.gzip.h5ad -r test/scRNAseq_downsample_0.gzip.h5ad -c cell_class1,seurat_clusters --continuous_columns_to_impute umapx,umapy -a 20000 -e 2 -m 30 -t 100 -n 8 
