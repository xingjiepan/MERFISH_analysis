import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_barcode_files_under_dir(barcode_dump_path):
    '''Load barcode files under a given directory.'''
    barcode_files = [os.path.join(barcode_dump_path, f) for f in os.listdir(barcode_dump_path)
                                                         if f.endswith('.h5')]
    df_list = [] 
    for bf in barcode_files:
        print(f'Load barcodes in {bf}')
        df_list.append(pd.read_hdf(bf))

    return pd.concat(df_list) 

def calc_pairwise_distances_between_points(P1:np.ndarray, P2:np.ndarray) -> np.ndarray:
    '''Calculate paiwise distances between two arrays of points.
    Return a matrix of pairwise distances.
    '''
    return np.sqrt(np.square(P1[:, np.newaxis, :] - P2[np.newaxis, :, :]).sum(axis=2))
    
def get_N_colocalized_barcodes(barcode_df, z1, z2, barcode_id, distance_threshold=2):
    '''Get the number of colocalized barcodes of a given barcode id between two z-planes.'''
    ids1 = np.logical_and(barcode_df['z'] == z1, barcode_df['barcode_id'] == barcode_id)
    ids2 = np.logical_and(barcode_df['z'] == z2, barcode_df['barcode_id'] == barcode_id)
    
    P1 = np.array(barcode_df[ids1][['x', 'y']])
    P2 = np.array(barcode_df[ids2][['x', 'y']])
    
    dist_mtx = calc_pairwise_distances_between_points(P1, P2)
    n1 = np.sum(np.sum(dist_mtx < distance_threshold, axis=1) > 0) # Number of rows
    n2 = np.sum(np.sum(dist_mtx < distance_threshold, axis=0) > 0) # Number of columns
    
    return n1, n2

def get_N_colocalized_barcodes_between_z_planes(barcode_df, Zs, distance_threshold):
    '''Get the numbers of coloalized barcodes between each pair
    of consecutive Z-planes.
    Return:
        n_barcodes: an 1D array of numbers of barcodes in each Z plane.
        n_colocalized_barcodes: an 1D array of numbers of colocalized barcodes between consecutive Z-planes.
    '''
    # Calculate the numbers of barcodes in each Z plane
    n_barcodes = []
    for z in Zs:
        n_barcodes.append(np.sum(barcode_df['z'] == z))
        
    # Calculate the numbers of colocalized barcodes between consecutive Z planes
    n_colocalized_barcodes = []
    barcode_ids = np.unique(barcode_df['barcode_id'])
    
    for i in range(len(Zs) - 1):
        z1 = Zs[i]
        z2 = Zs[i + 1]
        
        n_co = 0
    
        for bid in barcode_ids:
            n1, n2 = get_N_colocalized_barcodes(barcode_df, z1, z2, bid, 
                                                distance_threshold=distance_threshold)
            n_co += n1
                
        n_colocalized_barcodes.append(n_co)

    return np.array(n_barcodes), np.array(n_colocalized_barcodes)

def get_N_colocalized_barcodes_for_all_FOVs(barcode_dump_path, Zs, distance_threshold):
    '''Get the numbers of coloalized barcodes between each pair
    of consecutive Z-planes for all fields of views.
    Return:
        n_barcodes: a list of numbers of barcodes in each Z plane.
        n_colocalized_barcodes: a list of numbers of colocalized barcodes between consecutive Z-planes.
    '''
    barcode_files = [os.path.join(barcode_dump_path, f) for f in os.listdir(barcode_dump_path)
                                                         if f.endswith('.h5')]
    
    n_barcodes = np.zeros(len(Zs))
    n_colocalized_barcodes = np.zeros(len(Zs) - 1)
    
    for bf in barcode_files:
        print(f'Analyze barcodes in {bf}')
        barcode_df = pd.read_hdf(bf)
        
        nbs, ncbs = get_N_colocalized_barcodes_between_z_planes(barcode_df, Zs, 
                                                distance_threshold=distance_threshold)
        n_barcodes += nbs
        n_colocalized_barcodes += ncbs
        
    return n_barcodes, n_colocalized_barcodes

def plot_N_colocalized_barcodes_between_z_planes(Zs, n_barcodes, n_colocalized_barcodes):
    N_Zs = len(Zs)
    
    plt.bar(list(range(N_Zs - 1)), n_barcodes[:N_Zs - 1], label='total barcodes')
    plt.bar(list(range(N_Zs - 1)), n_colocalized_barcodes[:N_Zs - 1], label='colocalized barcodes')
    plt.legend()
    
    plt.xticks(list(range(N_Zs - 1)), Zs[:N_Zs - 1])
    plt.xlabel('Z')
    plt.ylabel('Number of barcodes')
    
    plt.show()

