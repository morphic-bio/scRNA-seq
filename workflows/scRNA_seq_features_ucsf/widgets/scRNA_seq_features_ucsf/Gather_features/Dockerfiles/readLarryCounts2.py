import argparse
import scanpy as sc
import os
import sys
import pyroe
import anndata as ad
import pandas as pd
import scipy.sparse as sp
import numpy as np
import h5py
from scipy.sparse import csr_matrix, lil_matrix

def add_cellbender_counts(sc_adata, cb_dir):
    counts_h5_path = os.path.join(cb_dir, 'cellbender_counts.h5')
    if not os.path.exists(counts_h5_path):
        raise FileNotFoundError(f"Required CellBender file not found in {cb_dir}")
    
    with h5py.File(counts_h5_path, 'r') as h5_file:
        cb_counts = h5_file['matrix/data'][:]
        cb_indices = h5_file['matrix/indices'][:]
        cb_indptr = h5_file['matrix/indptr'][:]
        cb_shape = h5_file['matrix/shape'][:]
        cb_barcodes = h5_file['matrix/barcodes'][:].astype(str) if 'matrix/barcodes' in h5_file else [f'cell_{i}' for i in range(len(cb_indptr) - 1)]
        
        cb_counts_matrix = csr_matrix((cb_counts, cb_indices, cb_indptr), shape=(len(cb_barcodes), cb_shape[1]))
    
    cb_adata = ad.AnnData(X=cb_counts_matrix, obs=pd.DataFrame(index=cb_barcodes))
    
    # Step 1: Align obs_names (barcodes)
    sc_adata.obs_names = sc_adata.obs_names.str.strip().astype(str)
    cb_adata.obs_names = cb_adata.obs_names.str.strip().astype(str)
    sc_adata.obs_names_make_unique()
    cb_adata.obs_names_make_unique()
    
    # Step 2: Align var_names (genes)
    sc_adata.var_names = sc_adata.var_names.str.strip().astype(str)
    cb_adata.var_names = cb_adata.var_names.str.strip().astype(str)
    sc_adata.var_names_make_unique()
    cb_adata.var_names_make_unique()
    
    sc_obs_indices = sc_adata.obs_names.get_indexer_for(cb_adata.obs_names)
    sc_var_indices = sc_adata.var_names.get_indexer_for(cb_adata.var_names)
    
    # Step 3: Create a sparse layer for denoised data
    denoised_layer = lil_matrix(sc_adata.shape, dtype=np.float32)
    
    # Step 4: Copy only matching barcodes and genes
    for i, row_idx in enumerate(sc_obs_indices):
        if row_idx == -1:
            continue  # Skip if barcode does not exist in sc_adata
        matching_genes = sc_var_indices[sc_var_indices != -1]
        denoised_layer[row_idx, matching_genes] = cb_counts_matrix[i, matching_genes].toarray().flatten()
    
    denoised_layer = denoised_layer.tocsr()
    sc_adata.layers['denoised'] = denoised_layer
    print(f"Added sparse 'denoised' layer with shape {denoised_layer.shape} to sc_adata")
    
    return sc_adata

def merge_feature_sc_data(feature_adata, sc_adata):
    """
    Merge two AnnData objects by intersecting obs_names and var_names, 
    and transfer sc_adata values to a layer in feature_adata.

    Parameters:
    -----------
    feature_adata : AnnData
        The feature-level AnnData object (e.g., targeted or multiomic dataset).
    sc_adata : AnnData
        The single-cell AnnData object to be merged with feature_adata.

    Returns:
    --------
    AnnData
        The updated feature AnnData object, with sc_counts layer added 
        and obs fields from sc_adata transferred to feature_adata.
    """
    # Step 1: Ensure obs_names and var_names are unique and clean
    feature_adata.obs_names = feature_adata.obs_names.str.strip().astype(str)
    sc_adata.obs_names = sc_adata.obs_names.str.strip().astype(str)
    
    feature_adata.obs_names_make_unique()
    sc_adata.obs_names_make_unique()
    
    feature_adata.var_names = feature_adata.var_names.str.strip().astype(str)
    sc_adata.var_names = sc_adata.var_names.str.strip().astype(str)
    
    feature_adata.var_names_make_unique()
    sc_adata.var_names_make_unique()

    # Step 2: Find intersection of obs_names (barcodes) and var_names (genes)
    common_obs_names = feature_adata.obs_names.intersection(sc_adata.obs_names)
    common_var_names = feature_adata.var_names.intersection(sc_adata.var_names)
    
    print(f"Number of common obs_names: {len(common_obs_names)}")
    print(f"Number of common var_names: {len(common_var_names)}")
    
    # Step 3: Subset feature_adata and sc_adata to these common indices
    feature_adata_subset = feature_adata[common_obs_names, common_var_names].copy()
    sc_adata_subset = sc_adata[common_obs_names, common_var_names].copy()
    
    # Step 4: Transfer sc_adata values to a new layer in feature_adata
    if isinstance(sc_adata_subset.X, np.ndarray):
        feature_adata_subset.layers['sc_counts'] = sc_adata_subset.X
    else:
        feature_adata_subset.layers['sc_counts'] = sc_adata_subset.X.toarray()
    
    # Step 5: Transfer all other obs fields from sc_adata to feature_adata
    for obs_column in sc_adata_subset.obs.columns:
        if obs_column not in feature_adata_subset.obs.columns:
            feature_adata_subset.obs[obs_column] = sc_adata_subset.obs[obs_column]
        else:
            print(f"Warning: Column '{obs_column}' already exists in feature_adata. Skipping merge for this column.")
    
    return feature_adata_subset
def merge_feature_sc_data(feature_adata, sc_adata):
    """
    Merge two AnnData objects by intersecting obs_names and var_names, 
    and transfer sc_adata values to a layer in feature_adata.

    Parameters:
    -----------
    feature_adata : AnnData
        The feature-level AnnData object (e.g., targeted or multiomic dataset).
    sc_adata : AnnData
        The single-cell AnnData object to be merged with feature_adata.

    Returns:
    --------
    AnnData
        The updated feature AnnData object, with sc_counts layer added 
        and obs fields from sc_adata transferred to feature_adata.
    """
    # Step 1: Ensure obs_names and var_names are unique and clean
    feature_adata.obs_names = feature_adata.obs_names.str.strip().astype(str)
    sc_adata.obs_names = sc_adata.obs_names.str.strip().astype(str)
    
    feature_adata.obs_names_make_unique()
    sc_adata.obs_names_make_unique()
    
    feature_adata.var_names = feature_adata.var_names.str.strip().astype(str)
    sc_adata.var_names = sc_adata.var_names.str.strip().astype(str)
    
    feature_adata.var_names_make_unique()
    sc_adata.var_names_make_unique()

    # Step 2: Find intersection of obs_names (barcodes) and var_names (genes)
    common_obs_names = feature_adata.obs_names.intersection(sc_adata.obs_names)
    common_var_names = feature_adata.var_names.intersection(sc_adata.var_names)
    
    print(f"Number of common obs_names: {len(common_obs_names)}")
    print(f"Number of common var_names: {len(common_var_names)}")
    
    # Step 3: Subset feature_adata and sc_adata to these common indices
    feature_adata_subset = feature_adata[common_obs_names, common_var_names].copy()
    sc_adata_subset = sc_adata[common_obs_names, common_var_names].copy()
    
    # Step 4: Transfer sc_adata values to a new layer in feature_adata
    if isinstance(sc_adata_subset.X, np.ndarray):
        feature_adata_subset.layers['sc_counts'] = sc_adata_subset.X
    else:
        feature_adata_subset.layers['sc_counts'] = sc_adata_subset.X.toarray()
    
    # Step 5: Transfer all other obs fields from sc_adata to feature_adata
    for obs_column in sc_adata_subset.obs.columns:
        if obs_column not in feature_adata_subset.obs.columns:
            feature_adata_subset.obs[obs_column] = sc_adata_subset.obs[obs_column]
        else:
            print(f"Warning: Column '{obs_column}' already exists in feature_adata. Skipping merge for this column.")
    
    return feature_adata_subset
def merge_adata_as_layers(adata1, adata2, layer1_name='', layer2_name='layer2'):
    """
    Merge two AnnData objects as layers with the obs_names being the union of the two.

    Parameters:
    -----------
    adata1 : AnnData
        The first AnnData object to merge.
    adata2 : AnnData
        The second AnnData object to merge.
    layer1_name : str
        Name of the layer corresponding to adata1.
    layer2_name : str
        Name of the layer corresponding to adata2.

    Returns:
    --------
    AnnData
        A new AnnData object containing both adata1 and adata2 as layers,
        with obs_names being the union of the two.
    """
    # Step 1: Ensure obs_names are clean and unique
    adata1.obs_names = adata1.obs_names.str.strip().astype(str)
    adata2.obs_names = adata2.obs_names.str.strip().astype(str)
    
    adata1.obs_names_make_unique()
    adata2.obs_names_make_unique()
    
    # Step 2: Create the union of obs_names
    union_obs_names = pd.Index(adata1.obs_names).union(pd.Index(adata2.obs_names))
    
    # Step 3: Create missing entries for adata1
    missing_obs_adata1 = list(set(union_obs_names) - set(adata1.obs_names))
    if len(missing_obs_adata1) > 0:
        # Create an empty AnnData with the missing observations (filled with zeros)
        empty_adata1 = ad.AnnData(X=np.zeros((len(missing_obs_adata1), adata1.shape[1])),
                                  obs=pd.DataFrame(index=missing_obs_adata1),
                                  var=adata1.var.copy())
        # Concatenate original adata1 with the new empty rows
        adata1 = ad.concat([adata1, empty_adata1], axis=0)

    # Step 4: Create missing entries for adata2
    missing_obs_adata2 = list(set(union_obs_names) - set(adata2.obs_names))
    if len(missing_obs_adata2) > 0:
        # Create an empty AnnData with the missing observations (filled with zeros)
        empty_adata2 = ad.AnnData(X=np.zeros((len(missing_obs_adata2), adata2.shape[1])),
                                  obs=pd.DataFrame(index=missing_obs_adata2),
                                  var=adata2.var.copy())
        # Concatenate original adata2 with the new empty rows
        adata2 = ad.concat([adata2, empty_adata2], axis=0)
    
    # Step 5: Reorder the obs to match the union_obs_names
    adata1 = adata1[union_obs_names, :]
    adata2 = adata2[union_obs_names, :]

    # Step 6: Create a new AnnData object with union_obs_names and the same var_names as adata1
    num_obs = len(union_obs_names)
    num_vars = len(adata1.var_names)
    merged_adata = ad.AnnData(X=None, 
                              obs=pd.DataFrame(index=union_obs_names), 
                              var=adata1.var.copy())
    
    # Initialize two layers with zeros (sparse or dense)
    if (layer1_name == ''):
        merged_adata.X = np.zeros((num_obs, num_vars)) 
    else: 
        merged_adata.layers[layer1_name] = np.zeros((num_obs, num_vars))
    merged_adata.layers[layer2_name] = np.zeros((num_obs, num_vars))

    # Step 7: Copy adata1's matrix into merged_adata's first layer
    if isinstance(adata1.X, np.ndarray):
        if (layer1_name == ''):
            merged_adata.X[:, :] = adata1.X
        else:
            merged_adata.layers[layer1_name][:, :] = adata1.X
    else:
        if (layer1_name == ''):
            merged_adata.X[:, :] = adata1.X.toarray()
        else:
            merged_adata.layers[layer1_name][:, :] = adata1.X.toarray()

    # Step 8: Copy adata2's matrix into merged_adata's second layer
    if isinstance(adata2.X, np.ndarray):
        merged_adata.layers[layer2_name][:, :] = adata2.X
    else:
        merged_adata.layers[layer2_name][:, :] = adata2.X.toarray()
    return merged_adata

def read_larry(directory):
    total_file=os.path.join(directory, 'features_matrix.mtx')
    bpath=os.path.join(directory, 'barcodes.txt')
    gpath=os.path.join(directory, 'features.txt')
    #read the matrices in
    print(total_file)
    x = sc.read_mtx(total_file)
    
    bc = pd.read_csv(bpath, names=['barcodes'])
    g = pd.read_csv(gpath, sep=' ', header=0, index_col=0)
    gene_names = g.index.tolist()
    barcode_names = bc.barcodes.tolist()

    bc = bc.set_index('barcodes')
    g = g.set_index(g.index)
    x=x.T
    x.var = g
    x.obs = bc
    adata = ad.AnnData(X=x.X)
    adata.var_names = gene_names
    adata.obs_names = barcode_names
    print(adata)
    return adata
# Initialize the parser
parser = argparse.ArgumentParser(description='Need to know where the base directory is for larry counts and the base directory for sc counts')

# Adding the arguments

parser.add_argument('--targeted-dir', required=True, type=str,
                    help='Specify the larry directory.')
parser.add_argument('--expression-dir', required=True, type=str,
                    help='Specify the expression directory.')                    
parser.add_argument('--sc-dir', required=True, type=str,
                    help='Specify the sc directory.')
parser.add_argument('--overwrite', action='store_true', help='Overwrite existing files')    

# Parse the arguments
args = parser.parse_args()
                    
overwrite = False                  
larry_base_directory = args.targeted_dir
sc_base_directory = args.sc_dir
expression_base_dir = args.expression_dir
overwrite = args.overwrite
#find all directries in sc_directory
print(f"Reading larry counts from {larry_base_directory}")
print(f"Reading expression counts from {expression_base_dir}")
print(f"Reading sc counts from {sc_base_directory}")
sub_directories = [d for d in os.listdir(sc_base_directory) if os.path.isdir(os.path.join(sc_base_directory, d))]
for dir in sub_directories:
    print(f"Reading counts from {dir}")
    larry_directory=os.path.join(larry_base_directory,dir)
    expression_directory=os.path.join(expression_base_dir,dir)
    sc_directory=os.path.join(sc_base_directory,dir)
    print(f"Reading larry counts from {larry_directory}")
    if os.path.exists(os.path.join(larry_directory, 'merged_features.h5ad')) and not overwrite:
        print(f"Skipping reading counts from  {larry_directory} as the output file already exists.")
        merged_adata = ad.read_h5ad(os.path.join(larry_directory, 'merged_features.h5ad'))
    else:
        larry_adata=read_larry(larry_directory)
        print(f"Reading expression counts from {expression_directory}")
        expression_adata=read_larry(expression_directory)
        merged_adata=merge_adata_as_layers(larry_adata, expression_adata, layer2_name='expression')
        merged_adata.write(os.path.join(larry_directory, 'merged_features.h5ad'))
        print(f"Saved merged AnnData object to {os.path.join(larry_directory, 'merged_features.h5ad')}")
    #print the shape of the merged AnnData object   
    print(merged_adata)
    methods_directories = [d for d in os.listdir(sc_directory) if os.path.isdir(os.path.join(sc_directory, d))]
    output_paths = []
    input_paths = []
    for method in methods_directories:
        if (method == 'salmon'):
            continue
        if (method == 'salmon' or method == 'piscem'):
            repos= [d for d in os.listdir(os.path.join(sc_directory, method)) if os.path.isdir(os.path.join(sc_directory, method, d))]
            for repo in repos:
                if (repo == 'splicei'):
                    output_paths.append(os.path.join(sc_directory,method, repo, 'merged_features.h5ad'))
                    input_paths.append(os.path.join(sc_directory, method, repo, 'unfiltered_counts.h5ad'))
        else:
            output_paths.append(os.path.join(sc_directory, method, 'merged_features.h5ad'))
            input_paths.append(os.path.join(sc_directory, method, 'unfiltered_counts.h5ad'))   
    for i in range(len(input_paths)):
        #check if input file exists
        input_path = input_paths[i]
        output_path = output_paths[i]
        if not os.path.exists(input_path):
            print(f"Skipping {input_path} as it does not exist")
            continue
        if os.path.exists(output_path) and not overwrite:
            print(f"Skipping reading counts from  {input_path} as the output file {output_path} already exists.")
            continue
        print(f"Reading counts from {input_path}")
        sc_adata=ad.read_h5ad(input_path)
        sc_cb_adata = add_cellbender_counts(sc_adata, os.path.join(os.path.dirname(input_path), 'cellbender'))
        print(sc_cb_adata) 
        cb_output_path = os.path.join(os.path.dirname(input_path), 'unfiltered_cb_counts.h5ad')
        sc_cb_adata.write(cb_output_path)
        sc_merged_adata = merge_feature_sc_data(merged_adata, sc_cb_adata)
        print(sc_merged_adata)
        merged_adata.write(output_path)
exit(0)


directories = [d for d in os.listdir(sc_directory) if os.path.isdir(os.path.join(sc_directory, d))]
for dir in directories:
    larry_directory=os.path.join(larry_base_directory,dir)
    print(f"Reading larry counts from {larry_directory}")
    larry_adata=read_larry(larry_directory)
    subdirs = [d for d in os.listdir(os.path.join(sc_directory, dir)) if os.path.isdir(os.path.join(sc_directory, dir, d))]
    for d in subdirs:
        #check if output file exists
        fullpath = os.path.join(sc_directory, dir, d, 'raw_features.h5ad')
        output_path = os.path.join(sc_directory,dir,d, 'merged_features.h5ad')
        # Skip if the fullpath does not exist
        if not os.path.exists(fullpath):
            print(f"Skipping {fullpath} as it does not exist")
            continue
        if os.path.exists(output_path):
            print(f"Skipping reading counts from  {fullpath} as the output file {output_path} already exists.")
            continue
        print(f"Reading counts from {fullpath}")
        sc_adata=ad.read_h5ad(fullpath)
        merged_adata = merge_anndata_sparse(larry_adata, sc_adata) 
        print(merged_adata)
        output_path = os.path.join(sc_directory,dir,d, 'merged_features.h5ad')
        merged_adata.write(output_path)
        print(f"Saved merged AnnData object to {output_path}")  