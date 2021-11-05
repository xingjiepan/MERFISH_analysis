#!/usr/bin/env python3
'''Generate the input files for hierarchical integration.
'''

from hierarchical_integration_lib import prepare_integration_inputs



if __name__ == '__main__':
    
    output_path = 'test/output'
    reference_adata_file = 'test/scRNAseq_downsample_0.gzip.h5ad'
    query_adata_file = 'test/merfish_downsample_0.gzip.h5ad'
    reference_cell_type_column = 'seurat_clusters'
    query_cell_type_column = 'de_novo_cluster'
    approximate_subset_size = 1000
    n_repeat_query = 2
    min_N_cells_per_cluster = 10

    prepare_integration_inputs(output_path, reference_adata_file, query_adata_file, 
        reference_cell_type_column, query_cell_type_column, approximate_subset_size,
        n_repeat_query, min_N_cells_per_cluster)
