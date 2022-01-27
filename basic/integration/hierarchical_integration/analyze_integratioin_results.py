#!/usr/bin/env python3
'''Analyze integration results.
'''

import os
from optparse import OptionParser
from collections import Counter

import numpy as np
import pandas as pd

import anndata



def load_integrated_cell_types(cell_type_files, mixing_score_files):
    '''Load the integrated cell types into a dictionary.
    The keys are cell ids and the values are the predicted cell types
    and predicted probabilities.
    '''
    ic_dict = {}
    
    for i, ct_file in enumerate(cell_type_files):
        print(f'Load integration result from {i}th file {ct_file}')
        
        ct_df = pd.read_csv(ct_file)
        ms_df = pd.read_csv(mixing_score_files[i])
        ct_df = ct_df.merge(ms_df, left_on='Unnamed: 0', right_on='Unnamed: 0', how='inner')

        cell_ids = [str(x) for x in ct_df.iloc[:,0]]
        p_types = list(ct_df['predicted.id'])
        p_probas = list(ct_df['prediction.score.max'])
        m_scores = list(ct_df['mixing_score'])

        for j in range(ct_df.shape[0]):
        
            cid = cell_ids[j]
            if not (cid in ic_dict):
                ic_dict[cid] = [[], [], []]
                
            ic_dict[cid][0].append(p_types[j]) # predicted cell type
            ic_dict[cid][1].append(p_probas[j]) # probability of the predicted type
            ic_dict[cid][2].append(m_scores[j]) # mixing score of the cell

    return ic_dict

def calculate_predicted_type_and_confidence(ic_dict):
    '''Calculate the confidence scores of each predicted cell.
    Return the results in a data frame.
    '''
    cell_ids = list(ic_dict.keys())
    predicted_types = []
    prediction_confidences = []
    mixing_scores = []
    
    for i, cid in enumerate(cell_ids):        
        cts = np.array(ic_dict[cid][0])
        cps = np.array(ic_dict[cid][1])
        
        counts = Counter(cts)
        mc_t = counts.most_common(1)[0][0]
        
        confidence = np.sum((cts == mc_t) * cps) / cps.shape[0]
        
        predicted_types.append(mc_t)
        prediction_confidences.append(confidence)
        mixing_scores.append(np.mean(ic_dict[cid][2]))
        
        if i % 10000 == 0:
            print(f'Processed {i} / {len(cell_ids)} cells.', end='\r')
    
    return pd.DataFrame({'id':cell_ids, 'predicted_type':predicted_types, 
                         'prediction_confidence':prediction_confidences,
                         'mixing_score':mixing_scores}).set_index('id')

def get_predicted_types_and_confidences(integration_path):
    '''Get the predicted types, confidences and mixing scores.'''
    # Find all the cell type files
    cell_type_files = []
    mixing_score_files = []
    for f in os.listdir(integration_path):
        subsets_path = os.path.join(integration_path, f, 'subsets')
        if os.path.exists(subsets_path):
            for ff in os.listdir(subsets_path):
                
                ct_file = os.path.join(subsets_path, ff, 'integrated', 'predicted_cell_types.csv')
                mixing_score_file = os.path.join(subsets_path, ff, 'integrated', 'query_cell_mixing_scores.csv')
                if os.path.exists(ct_file) and os.path.exists(mixing_score_file):
                    cell_type_files.append(ct_file)
                    mixing_score_files.append(mixing_score_file)

    # Load and calculate the confidences
    ic_dict = load_integrated_cell_types(cell_type_files, mixing_score_files)

    return calculate_predicted_type_and_confidence(ic_dict)
    
def analyze_integration_result(integration_path, adata_file_before_integration, 
        prediction_cell_type_col, mixing_threshold=100):
    '''Analyze the integration result and save it as an anndata file.'''
    
    # Load the integrated cell types and confidences
    ct_prediction_df = get_predicted_types_and_confidences(integration_path)
    ct_prediction_df.rename(columns={'predicted_type':prediction_cell_type_col, 
                                     'prediction_confidence':prediction_cell_type_col + '_confidence',
                                     'mixing_score':prediction_cell_type_col + '_mixing_score'}, inplace=True)

    # Set the cell types of low confidence cells to be pd.NA 
    col_filtered = np.array(ct_prediction_df[prediction_cell_type_col], dtype=str)
    col_filtered[ct_prediction_df[prediction_cell_type_col + '_mixing_score'] > mixing_threshold] = pd.NA
    ct_prediction_df[prediction_cell_type_col + '_filtered'] = col_filtered
    
    # Load the existing adata and add the integration information
    adata = anndata.read_h5ad(adata_file_before_integration)
    adata.obs = adata.obs.merge(ct_prediction_df, left_index=True, right_index=True, how='left')

    # Save the new adata file
    output_file = os.path.join(integration_path, 'integrated.h5ad')
    adata.write(output_file)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-t', '--mixing_threshold', dest='mixing_threshold', action='store', 
            type='float', default=100, help='The mixing score threshold to filter the predicted cell types.')

    (options, args) = parser.parse_args()

    integration_path = args[0]
    adata_file_before_integration = args[1]
    prediction_cell_type_col = args[2]
    mixing_threshold = options.mixing_threshold

    analyze_integration_result(integration_path, adata_file_before_integration, prediction_cell_type_col,
                               mixing_threshold)
