#!/bin/bash

./initialize_integration_project.py -i . -p test/test_integration_project -q test/merfish_downsample_0.gzip.h5ad -r test/scRNAseq_downsample_0.gzip.h5ad -c cell_class1,seurat_clusters --continuous_columns_to_impute umapx,umapy,PC0,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29 -a 20000 -e 2 -m 30 -t 100 -n 8 
