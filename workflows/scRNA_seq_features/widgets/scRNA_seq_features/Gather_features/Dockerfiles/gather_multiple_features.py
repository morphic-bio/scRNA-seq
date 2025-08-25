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


def add_called_features(features_dir: str,
                        called_features_dir: str,
                         counts_adata,
                         obs_key: str = "sample_assignment"
                         ):
    """
    Add a categorical obs column to `counts_adata` with sample assignments plus
    special categories: missing, low_support, ambiguous, doublets, filtered.

    Expected files in `called_features_dir`:
      - assignments.txt/tsv/csv      (2 cols: cell_barcode, feature_barcode)
      - ambiguous.txt/tsv/csv            (1st col: cell_barcode)
      - doublets.txt/tsv/csv             (1st col: cell_barcode)
      - missing_cells.txt/tsv/csv        (1st col: cell_barcode)
      - unassignable.txt/tsv/csv         (1st col: cell_barcode)  -> maps to 'low_support'
    Expected files in `features_dir`:
        - barcodes.txt/tsv/csv
        - features.txt/tsv/csv
        - matrix.mtx
    Behavior:
      - Cells in assignments are labeled by their id from features_list_file
      - If a sample name is missing for a feature, falls back to the feature_barcode.
      - Cells listed in the special files get those categories.
      - Any cell in counts_adata.obs_names not matched above → 'filtered'

    Returns:
      counts_adata with one or more categorical obs columns: obs_key] as a pandas.Categorical.
    """
    def _read_any(path):
        return pd.read_csv(path, sep=None, engine="python", header=None, comment="#")

    def _first_existing(directory,file_possibilities):
        for file in file_possibilities:
            q = os.path.join(directory, file)
            if os.path.exists(q):
                return q
        return None

    # --- locate files
    p_assign   = _first_existing(called_features_dir,["assignments.txt","assignments.tsv","assignments.csv"])
    p_names    = _first_existing(features_dir,["features.txt", "features.tsv", "features.csv"])
    p_amb      = _first_existing(called_features_dir,["ambiguous.txt", "ambiguous.tsv", "ambiguous.csv"])
    p_dbl      = _first_existing(called_features_dir,["doublets.txt", "doublets.tsv", "doublets.csv"])
    p_missing  = _first_existing(called_features_dir,["missing_cells.txt", "missing_cells.tsv", "missing_cells.csv"])
    p_low      = _first_existing(called_features_dir,["unassignable.txt", "unassignable.tsv", "unassignable.csv"])

    if p_assign is None:
        raise FileNotFoundError("Could not find features.assignments.(tsv|txt|csv) in features_dir")

    # --- read names mapping: feature_barcode -> sample_name
    feature_to_sample = {}
    idx_to_feature = None
    if p_names is not None:
        names_df = _read_any(p_names)
        # Expect a single column text file with no header - just feature names/barcode
        # Use first column as feature barcodes
        feature_barcodes = names_df.iloc[:, 0].astype(str).tolist()

        # Create idx_to_feature mapping from names file (1-based indexing)
        idx_to_feature = dict((i+1, feat) for i, feat in enumerate(feature_barcodes))


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
                "Assignments use 1-based indices but neither features.names.(txt|tsv|csv) nor allowed_features.(tsv|txt|csv) was found to map them."
            )
        assign_df["feature_barcode"] = assign_df["assigned"].map(idx_to_feature)
    else:
        # Already a feature barcode
        assign_df["feature_barcode"] = assign_df["assigned"].astype(str)


    assign_df["category"] = assign_df["feature_barcode"]

    # Deduplicate in case of repeated barcodes in file (keep first occurrence)
    assign_df = assign_df.drop_duplicates(subset=["cell_barcode"], keep="first")

    # --- read special categories (all optional files)
    def _read_barcode_list(path):
        """
        Read a single-column text file (optionally *.gz*) and return the
        first whitespace-separated token of each non-comment line as a `set`.
        """
        import gzip

        if path is None or not os.path.exists(path):
            return set()

        # Choose correct opener for plain vs gzip files
        opener = gzip.open if path.endswith(".gz") else open

        barcodes = set()
        try:
            with opener(path, "rt") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    # Keep the first token only (handles extra columns safely)
                    barcodes.add(line.split()[0])
        except Exception:
            # Any read/parse problem → treat as empty
            return set()

        return barcodes

    amb_set  = _read_barcode_list(p_amb)
    dbl_set  = _read_barcode_list(p_dbl)
    miss_set = _read_barcode_list(p_missing)   # -> "missing"
    low_set  = _read_barcode_list(p_low)       # -> "low_support"

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

# Uncomment the following lines only when you really need to debug environment variables
# if os.getenv("DEBUG_ENV"):
#     for key, value in os.environ.items():
#         print(f"{key}: {value}")

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
sample_names= ['30_KO_DE_XM','30_KO_ES','30_KO_PP1','30_KO_PP2','30_KO_S5_1','30_KO_S5_2','30_KO_S6_1','30_KO_S6_2']

for sample_name in sample_names:
    alignments_dir = f'/mnt/pikachu/storage/MSK-output-2/Alignments/{sample_name}/star/'
    counts_adata = ad.read_h5ad(f'{alignments_dir}/unfiltered_counts.h5ad')
    #check if the denoised layer exists
    if 'denoised' in counts_adata.layers:
        print("Denoised layer found")
    else:
        print("No denoised layer found - will add cellbender counts")
        cellbender_counts_h5ad = f'{alignments_dir}/cellbender/denoised_counts.h5ad'
        cellbender_counts_adata = ad.read_h5ad(cellbender_counts_h5ad)
        counts_adata = add_counts_layer_to_counts_adata(counts_adata, cellbender_counts_adata, 'denoised', 'denoised')
    features_dirs= [f'/storage/gene_features/{sample_name}',f'/storage/larry_features/{sample_name}']

    called_features_dirs= [f'/mnt/pikachu/storage/MSK-output-2/Alignments/{sample_name}/star/gene_features_em',f'/mnt/pikachu/storage/MSK-output-2/Alignments/{sample_name}/star/larry_features_em']
    anndata_feature_names=['gene_barcode','larry_barcode']
    for i,called_features_dir in enumerate(called_features_dirs):
        print(f"outputing features from {features_dirs[i]}")
        mtx_file = os.path.join(features_dirs[i], 'matrix.mtx')
        barcodes_file = os.path.join(features_dirs[i], 'barcodes.txt')
        features_file = os.path.join(features_dirs[i], 'features.txt')
        output_h5ad = os.path.join(features_dirs[i], 'features.h5ad')
        read_mtx_to_h5ad(barcodes_file, features_file, mtx_file, output_h5ad)
        counts_adata = add_called_features(features_dirs[i],    called_features_dir, counts_adata,anndata_feature_names[i])

    #mkdir if not exists
    os.makedirs(f'/mnt/pikachu/storage/MSK-output-2/Features/{sample_name}', exist_ok=True)    
    counts_adata.write(f'/mnt/pikachu/storage/MSK-output-2/Features/{sample_name}/features.h5ad')



