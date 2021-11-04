#!/usr/bin/env python3

import h5py
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

def extract_cell_boundaries_from_one_file(hdf5_feature_file):
    '''Return a dictionary. The keys are cell ids and the
    values are lists of boundaries at each z level.
    '''
    dict_cell_boundaries = {}
    
    with h5py.File(hdf5_feature_file, 'r') as f:
        cell_ids = list(f['featuredata'].keys())
        
        # Iterate through all cells
        for cid in cell_ids:
            z_ids = list(f['featuredata'][cid].keys())
            cell_boundaries = []
            
            # Iterate through all z-levels
            for zid in z_ids:

                if zid.startswith('zIndex_'):
                    z_boundaries = []
                    
                    if type(f['featuredata'][cid][zid]) is h5py._hl.group.Group:
                        
                        # Iterate through all polygons in the z-level
                        for polygon in f['featuredata'][cid][zid].keys():
                            z_boundaries.append(f['featuredata'][cid][zid][polygon]['coordinates'][0])
                    
                    cell_boundaries.append(z_boundaries)
                
            dict_cell_boundaries[cid] = cell_boundaries

    return dict_cell_boundaries

def extact_cell_boundaries_from_multiple_files(hdf5_feature_files):
    dict_cell_boundaries = {}
    for hdf5_feature_file in hdf5_feature_files:
        print(f'Extract cell boundaries from file {hdf5_feature_file}')
        dict_cell_boundaries_one_file = extract_cell_boundaries_from_one_file(hdf5_feature_file)
        
        dict_cell_boundaries = {**dict_cell_boundaries, **dict_cell_boundaries_one_file}
        
    return dict_cell_boundaries

def plot_polygons_at_z_level(dict_cell_boundaries, z):
    fig, ax = plt.subplots()
    
    for cell_id in dict_cell_boundaries:
        for polygon in dict_cell_boundaries[cell_id][z]:
            ax.plot(polygon[:,0], polygon[:,1])
            
    ax.set_aspect('equal')
    ax.grid(False)
    plt.show()
    
