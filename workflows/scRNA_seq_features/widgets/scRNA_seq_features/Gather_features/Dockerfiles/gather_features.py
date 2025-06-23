#! /usr/bin/env python3
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
import ast
import io
import scipy.io

def add_counts_layer_to_counts_adata(counts_adata, new_adata, source_layer_name, target_layer_name):
    """
    Adds a new sparse layer to counts_adata based on a specified layer from new_adata.
    
    This function is optimized for memory efficiency by avoiding dense matrix conversions.
    - For barcodes and genes present in both anndata objects, the values from the source layer are copied.
    - For barcodes or genes not present in the source, the values in the new layer will be 0 (implicit in sparse format).
    - A new boolean column `is_<target_layer_name>` is added to `counts_adata.obs` to indicate which
      barcodes were found in `new_adata`.

    Parameters
    ----------
    counts_adata : anndata.AnnData
        The target AnnData object to which the layer will be added.
    new_adata : anndata.AnnData
        The source AnnData object providing the data for the new layer.
    source_layer_name : str or None
        The name of the layer in `new_adata` to use as the source.
        If 'X' or None, `new_adata.X` is used.
    target_layer_name : str
        The name for the new layer in `counts_adata`.

    Returns
    -------
    anndata.AnnData
        The modified counts_adata object with the new sparse layer and obs column.
    """
    # Select the source data matrix
    if source_layer_name is None or source_layer_name == 'X':
        source_matrix = new_adata.X
        print(f"Using .X matrix from source anndata.")
    elif source_layer_name in new_adata.layers:
        source_matrix = new_adata.layers[source_layer_name]
        print(f"Using layer '{source_layer_name}' from source anndata.")
    else:
        raise ValueError(f"Source layer '{source_layer_name}' not found in source anndata.")

    # Find common observations (barcodes) and variables (genes)
    common_obs = counts_adata.obs_names.intersection(new_adata.obs_names)
    common_vars = counts_adata.var_names.intersection(new_adata.var_names)
    
    # Add a boolean obs column to indicate which barcodes were found in the source.
    obs_name = f"is_{target_layer_name}"
    counts_adata.obs[obs_name] = counts_adata.obs_names.isin(common_obs)
    
    print(f"Found {len(common_obs)} common barcodes and {len(common_vars)} common genes for layer '{target_layer_name}'.")
    print(f"Added boolean column '{obs_name}' to `counts_adata.obs`.")

    if len(common_obs) == 0 or len(common_vars) == 0:
        print(f"Warning: No common barcodes or genes found. Layer '{target_layer_name}' will be empty (all zeros).")
        counts_adata.layers[target_layer_name] = sp.csr_matrix(counts_adata.shape, dtype=source_matrix.dtype)
        return counts_adata

    # Get integer indices for mapping between matrices
    target_obs_indices = counts_adata.obs_names.get_indexer(common_obs)
    target_var_indices = counts_adata.var_names.get_indexer(common_vars)
    source_obs_indices = new_adata.obs_names.get_indexer(common_obs)
    source_var_indices = new_adata.var_names.get_indexer(common_vars)

    # Slice the source matrix to get data for common obs/vars, then convert to COO format
    source_data_common = source_matrix[source_obs_indices, :][:, source_var_indices]
    source_coo = source_data_common.tocoo(copy=False)

    # Remap the row and column indices from the source slice to the target matrix's coordinate system
    final_rows = target_obs_indices[source_coo.row]
    final_cols = target_var_indices[source_coo.col]

    # Build the new sparse layer directly from the remapped COO components
    new_layer_matrix = sp.coo_matrix((source_coo.data, (final_rows, final_cols)), shape=counts_adata.shape)

    # Add the new matrix as a layer to counts_adata, converting to CSR for efficiency.
    counts_adata.layers[target_layer_name] = new_layer_matrix.tocsr()

    return counts_adata

def filter_features(feature_adata, counts_adata):
    """
    Finds the top two features for each cell, their counts, and adds this
    information to counts_adata.obs.

    This function calculates the top and second-highest feature counts for each cell
    in feature_adata. It then adds the following columns to counts_adata.obs:
    - 'best_feature': The name of the feature with the highest count, if it is
      greater than the second-highest count. Otherwise, it's an empty string.
    - 'feature1_count': The count of the top feature.
    - 'feature2_count': The count of the second-highest feature.

    Parameters:
    -----------
    feature_adata : AnnData
        An AnnData object containing feature counts (cells x features).
    counts_adata : AnnData
        The AnnData object to which the new columns will be added.

    Returns:
    --------
    counts_adata : AnnData
        The updated counts_adata with the new columns in `obs`.
    """
    # Step 1: Find common barcodes
    common_barcodes = feature_adata.obs_names.intersection(counts_adata.obs_names)

    # Step 2: Select the correct data matrix. Prioritize the 'features' layer.
    if 'features' in feature_adata.layers:
        print("Using 'features' layer from feature_adata for filtering.")
        X_source = feature_adata.layers['features']
    else:
        print("Using .X from feature_adata for filtering.")
        X_source = feature_adata.X
    
    # Extract matrix
    X = X_source.toarray() if sp.issparse(X_source) else np.asarray(X_source)

    # Step 3: Find top two features and their counts
    if X.shape[1] > 1:
        # Get indices that would sort each row (cell) by counts in descending order
        sorted_indices = np.argsort(X, axis=1)[:, ::-1]
        
        # Indices of top and second-top features
        top_indices = sorted_indices[:, 0]
        second_top_indices = sorted_indices[:, 1]
        
        # Get counts of top and second-top features
        top_counts = X[np.arange(X.shape[0]), top_indices]
        second_top_counts = X[np.arange(X.shape[0]), second_top_indices]
        
        # Get names of the top features
        top_feature_names = feature_adata.var_names[top_indices]

    elif X.shape[1] == 1:
        top_counts = X[:, 0]
        second_top_counts = np.zeros_like(top_counts)
        top_feature_names = np.full(X.shape[0], feature_adata.var_names[0])
    else: # No features
        top_counts = np.zeros(X.shape[0])
        second_top_counts = np.zeros(X.shape[0])
        top_feature_names = np.full(X.shape[0], '')


    # Step 4: Create Series for new obs columns, indexed by feature_adata barcodes
    
    # Create 'best_feature' column
    # Condition: top count must be greater than second-top count
    best_feature_series = pd.Series(
        np.where(top_counts > second_top_counts, top_feature_names, ''),
        index=feature_adata.obs_names,
        name='best_feature'
    )
    
    # Create Series for top and second-top counts
    feature1_count_series = pd.Series(top_counts, index=feature_adata.obs_names, name='feature1_count')
    feature2_count_series = pd.Series(second_top_counts, index=feature_adata.obs_names, name='feature2_count')

    # Step 5: Add new columns to counts_adata.obs for common barcodes
    for obs_name, series_to_add in [('best_feature', best_feature_series),
                                     ('feature1_count', feature1_count_series),
                                     ('feature2_count', feature2_count_series)]:
        
        # Initialize column with appropriate empty value
        if pd.api.types.is_numeric_dtype(series_to_add):
            counts_adata.obs[obs_name] = np.nan
        else:
            counts_adata.obs[obs_name] = ''
            
        # Update the column for common barcodes with values from the series
        counts_adata.obs.loc[common_barcodes, obs_name] = series_to_add[common_barcodes]

    return counts_adata

def assign_feature(feature_adata, counts_adata, min_counts, obs_name):
    """
    Assigns the name of the feature (gene) with the highest count in each cell
    as a new column in the `obs` DataFrame of `feature_adata` and `counts_adata`.
    
    Parameters:
    -----------
    feature_adata : AnnData
        An AnnData object containing gene expression counts (cells x genes).
    counts_adata : AnnData
        Another AnnData object where the assigned feature will be copied.
    min_counts : float
        The minimum count required for the maximum value in a cell.
    obs_name : str
        The name of the new column to be added to `obs` for both AnnData objects.

    Returns:
    --------
    feature_adata : AnnData
        The updated feature_adata with the new column in `obs`.
    counts_adata : AnnData
        The updated counts_adata with the new column in `obs`.
    """
    # Step 1: Ensure consistent formatting of barcodes
    feature_adata.obs_names = feature_adata.obs_names.str.strip().astype(str)
    feature_adata.obs_names_make_unique()
    
    counts_adata.obs_names = counts_adata.obs_names.str.strip().astype(str)
    counts_adata.obs_names_make_unique()
    
    # Step 2: Extract the matrix from feature_adata
    X = feature_adata.X.toarray() if not isinstance(feature_adata.X, np.ndarray) else feature_adata.X
    
    # Get the maximum gene and second max value for each cell
    max_gene_indices = np.argmax(X, axis=1)  # Index of max value for each cell
    max_gene_values = np.max(X, axis=1)      # Value of max count for each cell
    X_temp = X.copy()
    X_temp[np.arange(X.shape[0]), max_gene_indices] = -np.inf
    second_max_gene_values = np.max(X_temp, axis=1)

    # Define which cells pass the feature assignment conditions
    condition_1 = max_gene_values >= min_counts
    condition_2 = second_max_gene_values <= 0.5 * max_gene_values
    valid_cells = condition_1 & condition_2

    # Assign the gene to the obs column for valid cells
    assigned_genes = [feature_adata.var_names[idx] if valid else '' for idx, valid in zip(max_gene_indices, valid_cells)]
    feature_adata.obs[obs_name] = assigned_genes

    print(f"Number of cells that met the condition: {np.sum(valid_cells)} / {feature_adata.shape[0]}")
    
    # Step 3: Align obs_names between feature_adata and counts_adata
    common_barcodes = feature_adata.obs_names.intersection(counts_adata.obs_names)
    missing_barcodes = feature_adata.obs_names.difference(counts_adata.obs_names)
    
    # Step 4: Align var_names (genes) to ensure proper column alignment
    common_var_names = feature_adata.var_names.intersection(counts_adata.var_names)
    
    if len(common_var_names) != len(counts_adata.var_names) or len(common_var_names) != len(feature_adata.var_names):
        feature_adata = feature_adata[:, common_var_names].copy()
        counts_adata = counts_adata[:, common_var_names].copy()
    
    # Step 5: Handle NaNs properly
    counts_adata.obs = counts_adata.obs.fillna('').astype(str)
    
    # Step 6: Add the assigned feature to counts_adata
    counts_adata.obs[obs_name] = feature_adata.obs[obs_name].astype(str)
    
    return feature_adata, counts_adata


def read_larry(directory):
    total_file=os.path.join(directory, 'features_matrix.mtx')
    bpath=os.path.join(directory, 'barcodes.txt')
    gpath=os.path.join(directory, 'features.txt')
    #read the matrices in
    print(total_file)
    
    with open(total_file, 'r') as f:
        lines = f.readlines()

    # Find the banner and the first non-comment line (which should be the dimension line)
    banner = lines[0]
    data_lines = [line for line in lines[1:] if not line.startswith('%')]
    
    # Reconstruct the file content in memory for mmread
    reconstructed_content = banner + ''.join(data_lines)
    
    matrix = scipy.io.mmread(io.StringIO(reconstructed_content))
    
    bc = pd.read_csv(bpath, names=['barcodes'])
    # We only need the index from the features file, not the data
    g = pd.read_csv(gpath, sep=' ', header=0, index_col=0)
    gene_names = g.index.tolist()
    barcode_names = bc.barcodes.tolist()

    # The matrix is (genes, cells), so we transpose it to (cells, genes).
    # Convert to CSR format for efficient storage and compatibility with h5ad writing.
    adata = ad.AnnData(X=matrix.T.tocsr())
    adata.var_names = gene_names
    adata.obs_names = barcode_names

    print(adata.shape)
    return adata


def merge_features_with_counts(features_adata, counts_adata):
    # Find intersection of barcodes that exist in both datasets
    common_barcodes = features_adata.obs_names.intersection(counts_adata.obs_names)

    # Add a boolean column to indicate if the barcode is in counts_adata
    features_adata.obs['is_expressed'] = features_adata.obs_names.isin(common_barcodes)
    counts_adata.obs['is_featured'] = counts_adata.obs_names.isin(common_barcodes)

    if len(common_barcodes) == 0:
        print("Warning: No common barcodes found between features_adata and counts_adata")
        return features_adata

    print(f"Found {len(common_barcodes)} common barcodes out of {len(features_adata.obs_names)} features barcodes")

    # Get the obs data from counts_adata for the common barcodes
    counts_obs_common = counts_adata.obs.loc[common_barcodes]

    # Iterate over columns in the counts metadata
    for col in counts_obs_common.columns:
        # Assign the data to features_adata.obs for common barcodes.
        # This will add the column if it doesn't exist (filling with NaN for non-common cells)
        # or update the existing values in the column for common cells.
        features_adata.obs.loc[common_barcodes, col] = counts_obs_common[col]

        # After assignment, if a column has an 'object' dtype, it may contain mixed
        # types (e.g., bools and NaNs), which causes issues when writing to h5ad.
        # Convert such columns to string to ensure they can be saved correctly.
        if features_adata.obs[col].dtype == 'object':
            features_adata.obs[col] = features_adata.obs[col].astype(str).fillna('')


    return features_adata, counts_adata



def merge_features(a1: ad.AnnData,
                   a2: ad.AnnData,
                   layer_names=('features', 'gex_features')) -> ad.AnnData:
    """
    Merge two AnnData objects with union of observations and variables.
    Each dataset goes into a separate layer.
    """
    # Get union of observations and variables
    obs_union = pd.Index(a1.obs_names).union(pd.Index(a2.obs_names))
    var_union = pd.Index(a1.var_names).union(pd.Index(a2.var_names))
    
    # Create new AnnData with union dimensions
    merged = ad.AnnData(
        X=sp.csr_matrix((len(obs_union), len(var_union))),
        obs=pd.DataFrame(index=obs_union),
        var=pd.DataFrame(index=var_union)
    )
    
    # Helper function to reindex sparse matrix
    def reindex_sparse_matrix(matrix, old_obs, old_var, new_obs, new_var):
        """Reindex sparse matrix to match new dimensions with union of obs/var"""
        # Get indexers for reindexing
        obs_indexer = new_obs.get_indexer(old_obs)
        var_indexer = new_var.get_indexer(old_var)
        
        # Create new sparse matrix with union dimensions
        new_matrix = sp.lil_matrix((len(new_obs), len(new_var)))
        
        # Fill in the data where indices exist
        valid_obs = obs_indexer >= 0
        valid_var = var_indexer >= 0
        
        if valid_obs.any() and valid_var.any():
            obs_idx_valid = obs_indexer[valid_obs]
            var_idx_valid = var_indexer[valid_var]
            
            # Extract submatrix for valid indices
            sub_matrix = matrix[valid_obs, :][:, valid_var]
            
            # Assign to new matrix
            new_matrix[np.ix_(obs_idx_valid, var_idx_valid)] = sub_matrix
        
        return new_matrix.tocsr()
    
    # Add data from a1 to first layer
    merged.layers[layer_names[0]] = reindex_sparse_matrix(
        a1.X, a1.obs_names, a1.var_names, obs_union, var_union
    )
    
    # Add data from a2 to second layer  
    merged.layers[layer_names[1]] = reindex_sparse_matrix(
        a2.X, a2.obs_names, a2.var_names, obs_union, var_union
    )
    
    # Merge obs and var metadata
    # For obs: keep all columns from both, fill missing with NaN
    obs_cols = set(a1.obs.columns) | set(a2.obs.columns)
    for col in obs_cols:
        merged.obs[col] = pd.Series(index=obs_union, dtype='object')
        if col in a1.obs.columns:
            merged.obs.loc[a1.obs_names, col] = a1.obs[col]
        if col in a2.obs.columns:
            merged.obs.loc[a2.obs_names, col] = a2.obs[col]
    
    # For var: keep all columns from both, fill missing with NaN  
    var_cols = set(a1.var.columns) | set(a2.var.columns)
    for col in var_cols:
        merged.var[col] = pd.Series(index=var_union, dtype='object')
        if col in a1.var.columns:
            merged.var.loc[a1.var_names, col] = a1.var[col]
        if col in a2.var.columns:
            merged.var.loc[a2.var_names, col] = a2.var[col]
    
    return merged

#print out all the env variables
for key, value in os.environ.items():
    print(f"{key}: {value}")

#find envs
features_dir_str = os.getenv('features_dirs')
alignments_dir = os.getenv('aligndir')
features_h5ad= os.getenv('features_h5ad')
merged_features_h5ad= os.getenv('merged_features_h5ad')
merged_counts_h5ad= os.getenv('merged_counts_h5ad')
input_counts_h5ad= os.getenv('input_counts_h5ad')
features_gex= os.getenv('features_gex_dirs')
cellbender_file= os.getenv('cellbender_file')
cb_original_layer= os.getenv('cb_original_layer')
cb_final_layer= os.getenv('cb_final_layer')
feature_assignment_str= os.getenv('feature_assignment')
overwrite = True
#overwrite= os.getenv('overwrite')
#if features_h5ad is not set set it to features.h5ad
if not features_h5ad:
    features_h5ad = 'features.h5ad'
    print(f"features_h5ad not set, setting it to {features_h5ad}")
if not merged_features_h5ad:
    merged_features_h5ad = 'merged_features.h5ad'
    print(f"merged_features_h5ad not set, setting it to {merged_features_h5ad}")
if not merged_counts_h5ad:
    merged_counts_h5ad = 'merged_counts.h5ad'
    print(f"merged_counts_h5ad not set, setting it to {merged_counts_h5ad}")
if not input_counts_h5ad:
    input_counts_h5ad = 'unfiltered_counts.h5ad'
    print(f"input_counts_h5ad not set, setting it to {input_counts_h5ad}")
features_dir = ast.literal_eval(features_dir_str) if features_dir_str else []

#get the features
for features_dir in features_dir:
    #get the sample which is the last part of the features_dir
    #use basename to get the sample
    sample = os.path.basename(features_dir)
    print(f"sample: {sample}")
    print(f"features_gex: {features_gex}")
    

    feature_adata = read_larry(features_dir)
    #print the shape of feature_adata
    print(f"shape of feature_adata: {feature_adata.shape}")
    #check if features_gex_sample_dir exists and is a directory
    if features_gex:
        print("features_gex is not empty")
        features_gex_sample_dir = os.path.join(features_gex, sample)
        if os.path.exists(features_gex_sample_dir) and os.path.isdir(features_gex_sample_dir):
            print ("reading features_gex_sample_dir: ", features_gex_sample_dir)
            features_gex_adata = read_larry(features_gex_sample_dir)
            print("merging features_gex_adata with feature_adata")
            feature_adata = merge_features(feature_adata, features_gex_adata)
    #write the feature_adata to a h5ad file
    #method is the subdirectory of the alignments_dir
    alignments_sample_dir = os.path.join(alignments_dir, sample)
    methods = [d for d in os.listdir(alignments_sample_dir) if os.path.isdir(os.path.join(alignments_sample_dir, d))]
    for method in methods:
        counts_h5ad = os.path.join(alignments_sample_dir, method, input_counts_h5ad)
        counts_adata = ad.read_h5ad(counts_h5ad)
        feature_adata, counts_adata = merge_features_with_counts(feature_adata, counts_adata)
        counts_adata=filter_features(feature_adata, counts_adata)
        feature_adata.write(os.path.join(alignments_sample_dir, method,merged_features_h5ad ))
        #add the cellbender counts to the counts_adata
        cellbender_counts_h5ad = os.path.join(alignments_sample_dir, method, cellbender_file)
        cellbender_counts_adata = ad.read_h5ad(cellbender_counts_h5ad)
        counts_adata = add_counts_layer_to_counts_adata(counts_adata, cellbender_counts_adata, cb_original_layer, cb_final_layer)
        counts_adata.write(os.path.join(alignments_sample_dir, method, merged_counts_h5ad))

