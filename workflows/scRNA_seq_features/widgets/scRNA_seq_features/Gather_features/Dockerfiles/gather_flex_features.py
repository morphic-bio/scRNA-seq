#! /usr/bin/env python3
from __future__ import annotations
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


def assign_flex_features(features_dir: str,
                         counts_adata,
                         obs_key: str = "sample_assignment"):
    """
    Add a categorical obs column to `counts_adata` with sample assignments plus
    special categories: missing, low_support, ambiguous, doublets, filtered.

    Expected files in `features_dir`:
      - features.assignments.tsv          (2 cols: cell_barcode, assigned_index [1-based])
      - allowed_features.tsv              (1 col: feature_barcode in the same order used for indices)
      - features.names.txt/tsv/csv        (2 cols: feature_barcode, sample_name)  # typo "nanes" handled
      - features.ambiguous.tsv            (1st col: cell_barcode)
      - features.doublets.tsv             (1st col: cell_barcode)
      - features.missing_cells.txt        (1st col: cell_barcode)
      - features.unassignable.txt         (1st col: cell_barcode)  -> maps to 'low_support'

    Behavior:
      - Cells in assignments are labeled by their *sample_name* from features.names.*.
        If a sample name is missing for a feature, falls back to the feature_barcode.
      - Cells listed in the special files get those categories.
      - Any cell in counts_adata.obs_names not matched above → 'filtered'.

    Returns:
      counts_adata with counts_adata.obs[obs_key] as a pandas.Categorical.
    """
    def _read_any(path):
        return pd.read_csv(path, sep=None, engine="python", header=None, comment="#")

    def _first_existing(*candidates):
        for p in candidates:
            q = os.path.join(features_dir, p)
            if os.path.exists(q):
                return q
        return None

    # --- locate files
    p_assign   = _first_existing("features.assignments.tsv", "features.assignments.txt", "features.assignments.csv")
    p_allowed  = _first_existing("allowed_features.tsv", "allowed_features.txt", "allowed_features.csv")
    p_names    = _first_existing("features.names.txt", "features.names.tsv", "features.names.csv",
                                 "features.nanes.txt", "features.nanes.tsv", "features.nanes.csv")  # typo tolerant
    p_amb      = _first_existing("features.ambiguous.tsv", "features.ambiguous.txt", "features.ambiguous.csv")
    p_dbl      = _first_existing("features.doublets.tsv", "features.doublets.txt", "features.doublets.csv")
    p_missing  = _first_existing("features.missing_cells.txt", "features.missing_cells.tsv", "features.missing_cells.csv")
    p_low      = _first_existing("features.unassignable.txt", "features.unassignable.tsv", "features.unassignable.csv")

    if p_assign is None:
        raise FileNotFoundError("Could not find features.assignments.(tsv|txt|csv) in features_dir")

    # --- read allowed feature list (index -> feature_barcode)
    idx_to_feature = None
    if p_allowed is not None:
        allowed_df = _read_any(p_allowed)
        allowed_df = allowed_df.iloc[:, :1]
        allowed_df.columns = ["feature_barcode"]
        # 1-based indexing in assignments → position in this list + 1
        idx_to_feature = dict((i+1, feat) for i, feat in enumerate(allowed_df["feature_barcode"].astype(str).tolist()))

    # --- read names mapping: feature_barcode -> sample_name
    feature_to_sample = {}
    if p_names is not None:
        names_df = _read_any(p_names)
        if names_df.shape[1] < 2:
            raise ValueError(f"{os.path.basename(p_names)} must have at least 2 columns: feature_barcode, sample_name")
        names_df = names_df.iloc[:, :2]
        names_df.columns = ["feature_barcode", "sample_name"]
        feature_to_sample = dict(zip(names_df["feature_barcode"].astype(str), names_df["sample_name"].astype(str)))

    # --- read assignments (cell_barcode, index_or_feature)
    assign_df = _read_any(p_assign)
    if assign_df.shape[1] < 2:
        raise ValueError(f"{os.path.basename(p_assign)} must have 2+ columns: cell_barcode, assigned_index")
    assign_df = assign_df.iloc[:, :2].copy()
    assign_df.columns = ["cell_barcode", "assigned"]
    assign_df["cell_barcode"] = assign_df["cell_barcode"].astype(str)

    # Is the second column numeric indices (1-based) or already feature barcodes?
    def _all_int_like(s):
        try:
            # allow whitespace, plus/minus; reject NaN
            return pd.to_numeric(s, errors="coerce").notna().all()
        except Exception:
            return False

    if _all_int_like(assign_df["assigned"]):
        assign_df["assigned"] = assign_df["assigned"].astype(int)
        if idx_to_feature is None:
            raise FileNotFoundError(
                "Assignments use 1-based indices but allowed_features.(tsv|txt|csv) was not found to map them."
            )
        assign_df["feature_barcode"] = assign_df["assigned"].map(idx_to_feature)
    else:
        # Already a feature barcode
        assign_df["feature_barcode"] = assign_df["assigned"].astype(str)

    # Map feature_barcode → sample_name (fallback to feature_barcode if name missing)
    if feature_to_sample:
        assign_df["category"] = assign_df["feature_barcode"].map(feature_to_sample).fillna(assign_df["feature_barcode"])
    else:
        # No names file: use feature_barcode as category
        assign_df["category"] = assign_df["feature_barcode"]

    # Deduplicate in case of repeated barcodes in file (keep first occurrence)
    assign_df = assign_df.drop_duplicates(subset=["cell_barcode"], keep="first")

    # --- read special categories (all optional files)
    def _read_barcode_list(path):
        if path is None: return set()
        df = _read_any(path)
        if df.shape[1] == 0:
            return set()
        return set(df.iloc[:, 0].astype(str).tolist())

    amb_set = _read_barcode_list(p_amb)
    dbl_set = _read_barcode_list(p_dbl)
    miss_set = _read_barcode_list(p_missing)     # -> "missing"
    low_set = _read_barcode_list(p_low)          # -> "low_support"

    # --- assemble final categorical series aligned to counts_adata.obs_names
    obs_index = counts_adata.obs_names.astype(str)
    result = pd.Series(index=obs_index, data=pd.NA, dtype="object")

    # Assign categories by priority:
    # 1) explicit special categories
    if len(miss_set):        result.loc[result.index.isin(miss_set)] = "missing"
    if len(low_set):         result.loc[result.index.isin(low_set)] = "low_support"
    if len(amb_set):         result.loc[result.index.isin(amb_set)] = "ambiguous"
    if len(dbl_set):         result.loc[result.index.isin(dbl_set)] = "doublets"

    # 2) sample assignments (do not overwrite special categories above)
    assigned_map = dict(zip(assign_df["cell_barcode"], assign_df["category"]))
    to_assign_mask = result.isna() & result.index.to_series().isin(assigned_map.keys())
    if to_assign_mask.any():
        result.loc[to_assign_mask] = result.index[to_assign_mask].map(assigned_map)

    # 3) everything else present in counts_adata → filtered
    result = result.fillna("filtered")

    # Build ordered category list: samples first (sorted), then special buckets
    sample_names = sorted(set(result.unique()) - {"missing", "low_support", "ambiguous", "doublets", "filtered"})
    categories = sample_names + ["missing", "low_support", "ambiguous", "doublets", "filtered"]
    result = pd.Categorical(result, categories=categories, ordered=False)

    # Attach to AnnData
    counts_adata.obs[obs_key] = result

    return counts_adata


def read_mtx_to_h5ad(barcodes_file: str, features_file: str, mtx_file: str, output_h5ad: str):
    """
    Read Matrix Market format files (barcodes, features, mtx) and write to h5ad format.
    
    This function reads the standard 10X Genomics output format:
    - barcodes.tsv(.gz): Cell barcodes (one per line)
    - features.tsv(.gz) or genes.tsv(.gz): Gene information (gene_id, gene_symbol, ...)
    - matrix.mtx(.gz): Sparse count matrix in Matrix Market format
    
    Parameters
    ----------
    barcodes_file : str
        Path to the barcodes file (can be .gz compressed)
    features_file : str
        Path to the features/genes file (can be .gz compressed)
    mtx_file : str
        Path to the matrix.mtx file (can be .gz compressed)
    output_h5ad : str
        Path for the output h5ad file
        
    Returns
    -------
    anndata.AnnData
        The loaded AnnData object
    """
    print(f"Reading matrix files:")
    print(f"  Barcodes: {barcodes_file}")
    print(f"  Features: {features_file}")
    print(f"  Matrix: {mtx_file}")
    
    # Check if files exist
    for file_path in [barcodes_file, features_file, mtx_file]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
    
    # Read barcodes
    print("Reading barcodes...")
    if barcodes_file.endswith('.gz'):
        import gzip
        with gzip.open(barcodes_file, 'rt') as f:
            barcodes = [line.strip() for line in f]
    else:
        with open(barcodes_file, 'r') as f:
            barcodes = [line.strip() for line in f]
    
    print(f"Found {len(barcodes)} barcodes")
    
    # Read features/genes
    print("Reading features...")
    try:
        # Try reading with pandas to handle different formats
        if features_file.endswith('.gz'):
            features_df = pd.read_csv(features_file, sep='\t', header=None, compression='gzip')
        else:
            features_df = pd.read_csv(features_file, sep='\t', header=None)
    except Exception as e:
        raise ValueError(f"Error reading features file {features_file}: {e}")
    
    # Handle different feature file formats
    if features_df.shape[1] >= 2:
        gene_ids = features_df.iloc[:, 0].astype(str).tolist()
        gene_symbols = features_df.iloc[:, 1].astype(str).tolist()
        print(f"Found {len(gene_ids)} genes with IDs and symbols")
    elif features_df.shape[1] == 1:
        gene_ids = features_df.iloc[:, 0].astype(str).tolist()
        gene_symbols = gene_ids.copy()  # Use IDs as symbols if only one column
        print(f"Found {len(gene_ids)} genes (using IDs as symbols)")
    else:
        raise ValueError(f"Features file must have at least 1 column, found {features_df.shape[1]}")
    
    # Read matrix
    print("Reading matrix...")
    try:
        matrix = scipy.io.mmread(mtx_file)
        # Convert to CSR format for efficiency
        matrix = matrix.tocsr()
        print(f"Matrix shape: {matrix.shape}")
    except Exception as e:
        raise ValueError(f"Error reading matrix file {mtx_file}: {e}")
    
    # Verify dimensions match
    if matrix.shape[0] != len(gene_ids):
        raise ValueError(f"Matrix has {matrix.shape[0]} rows but found {len(gene_ids)} genes")
    if matrix.shape[1] != len(barcodes):
        raise ValueError(f"Matrix has {matrix.shape[1]} columns but found {len(barcodes)} barcodes")
    
    # Create AnnData object
    print("Creating AnnData object...")
    
    # Create var DataFrame (genes/features)
    var_df = pd.DataFrame({
        'gene_ids': gene_ids,
        'feature_types': ['Gene Expression'] * len(gene_ids)  # Default feature type
    })
    var_df.index = gene_symbols
    
    # Create obs DataFrame (cells/barcodes)
    obs_df = pd.DataFrame(index=barcodes)
    
    # Create AnnData object (transpose matrix since mtx format is genes x cells)
    adata = ad.AnnData(X=matrix.T, obs=obs_df, var=var_df)
    
    # Add some basic information
    adata.var['gene_ids'] = gene_ids
    adata.var['gene_symbols'] = gene_symbols
    
    print(f"Created AnnData object with {adata.n_obs} cells and {adata.n_vars} genes")
    
    # Write to h5ad file
    print(f"Writing to {output_h5ad}...")
    adata.write(output_h5ad)
    print("Done!")
    
    return adata

#print out all the env variables
for key, value in os.environ.items():
    print(f"{key}: {value}")

#find envs
features_dir_str = os.getenv('features_dirs')
alignments_dir = os.getenv('aligndir')
features_h5ad= os.getenv('features_h5ad')
merged_features_h5ad= os.getenv('output_features_name')
merged_counts_h5ad= os.getenv('output_counts_name')
input_counts_h5ad= os.getenv('input_counts_name')
cellbender_file= os.getenv('cellbender_file')
cb_original_layer= os.getenv('cb_original_layer')
cb_final_layer= os.getenv('cb_final_layer')
overwrite = True

#for testing lets assign the above


features_dir ='/storage/scRNAseq_output/features/SC2300771'
barcodes_file = os.path.join(features_dir, 'barcodes.tsv')
features_file = os.path.join(features_dir, 'features.tsv')
mtx_file = os.path.join(features_dir, 'matrix.mtx')
output_h5ad = os.path.join(features_dir, 'features.h5ad')
read_mtx_to_h5ad(barcodes_file, features_file, mtx_file, output_h5ad)
exit()
counts_adata = ad.read_h5ad(input_counts_h5ad)
cellbender_counts_h5ad = os.path.join(alignments_sample_dir, method, cellbender_file)
cellbender_counts_adata = ad.read_h5ad(cellbender_counts_h5ad)
counts_adata = add_counts_layer_to_counts_adata(counts_adata, cellbender_counts_adata, cb_original_layer, cb_final_layer)
counts_adata = assign_flex_features(features_dir, counts_adata)
counts_adata.write(os.path.join(alignments_sample_dir, method, merged_counts_h5ad))


