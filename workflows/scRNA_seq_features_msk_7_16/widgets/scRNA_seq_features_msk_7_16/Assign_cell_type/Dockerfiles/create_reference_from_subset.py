import os
import anndata as ad
import pandas as pd
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import argparse

def create_reference_from_subset(celltype_column='celltype'):
    """
    Create a reference dataset by combining filtered adata files with celltypes.
    
    This function:
    1. Finds all merged_counts_with_celltypes.h5ad files created by add_msk_village.py
    2. Filters each by singlet_filtered=True and non-empty celltype
    3. Adds unique suffixes to barcodes to distinguish experiments
    4. Concatenates all filtered datasets
    5. Saves as sparse matrix to combined_reference.h5ad
    
    Parameters:
    -----------
    celltype_column : str
        Name of the obs column containing cell type information (default: 'celltype')
    """
    
    # Set up paths
    alignments_dir = os.getenv('alignments_dir')
    if not alignments_dir:
        alignments_dir = '/mnt/pikachu/storage/scRNAseq_output/Alignments'
    
    output_path = '/mnt/pikachu/scRNA-seq/data/combined_reference.h5ad'
    
    print(f"Using '{celltype_column}' column for cell type filtering")
    
    # Find all adata files with celltypes
    adata_files = []
    
    print("Searching for adata files with celltypes...")
    for sample in os.listdir(alignments_dir):
        sample_dir = os.path.join(alignments_dir, sample)
        if not os.path.isdir(sample_dir):
            continue
            
        for method in os.listdir(sample_dir):
            method_dir = os.path.join(sample_dir, method)
            if not os.path.isdir(method_dir):
                continue
                
            celltype_file = os.path.join(method_dir, 'merged_counts_with_celltypes.h5ad')
            if os.path.exists(celltype_file):
                adata_files.append({
                    'path': celltype_file,
                    'sample': sample,
                    'method': method,
                    'suffix': sample  # Changed from f"{sample}_{method}" to just sample
                })
                print(f"Found: {sample}/{method}")
    
    if not adata_files:
        print("No adata files with celltypes found!")
        return
    
    print(f"\nProcessing {len(adata_files)} files...")
    
    # Process and filter each file
    filtered_adatas = []
    
    for file_info in adata_files:
        print(f"\nProcessing {file_info['sample']}/{file_info['method']}...")
        
        # Load the data
        adata = ad.read_h5ad(file_info['path'])
        
        print(f"  Original shape: {adata.shape}")
        
        # Filter by singlet_filtered = True
        if 'singlet_filtered' in adata.obs.columns:
            singlet_mask = adata.obs['singlet_filtered'] == True
            adata_singlets = adata[singlet_mask].copy()
            print(f"  After singlet filtering: {adata_singlets.shape}")
        else:
            print("  Warning: 'singlet_filtered' column not found, skipping singlet filter")
            adata_singlets = adata.copy()
        
        # Filter by non-empty celltype, non-NaN celltype, and exclude 'nan' text
        if celltype_column in adata_singlets.obs.columns:
            celltype_mask = (adata_singlets.obs[celltype_column] != '') & \
                           (adata_singlets.obs[celltype_column].notna()) & \
                           (~adata_singlets.obs[celltype_column].isna()) & \
                           (adata_singlets.obs[celltype_column].str.lower() != 'nan')
            adata_filtered = adata_singlets[celltype_mask].copy()
            print(f"  After {celltype_column} filtering (excluding empty, NaN, and 'nan' text): {adata_filtered.shape}")
        else:
            print(f"  Warning: '{celltype_column}' column not found, skipping celltype filter")
            adata_filtered = adata_singlets.copy()
        
        if adata_filtered.shape[0] == 0:
            print("  No cells remaining after filtering, skipping...")
            continue
        
        # Apply lognorm tp10k normalization
        adata_filtered = apply_lognorm_tp10k(adata_filtered)
        
        # Add suffix to barcodes to make them unique
        suffix = file_info['suffix']
        new_obs_names = [f"{barcode}_{suffix}" for barcode in adata_filtered.obs_names]
        adata_filtered.obs_names = new_obs_names
        
        # Add sample and method information to obs
        adata_filtered.obs['sample'] = file_info['sample']
        adata_filtered.obs['method'] = file_info['method']
        
        # Ensure X matrix is CSR format (after normalization)
        if hasattr(adata_filtered.X, 'toarray'):
            if not sp.isspmatrix_csr(adata_filtered.X):
                adata_filtered.X = adata_filtered.X.tocsr()
        
        filtered_adatas.append(adata_filtered)
        print(f"  Added {adata_filtered.shape[0]} normalized cells to reference")
    
    if not filtered_adatas:
        print("No cells remaining after filtering any files!")
        return
    
    print(f"\nConcatenating {len(filtered_adatas)} filtered datasets...")
    
    # Concatenate all filtered datasets
    combined_adata = ad.concat(
        filtered_adatas,
        axis=0,  # concatenate along cell axis
        join='outer',  # include all genes
        merge='unique',  # handle duplicate gene names
        uns_merge='unique'
    )
    
    # Ensure the final matrix is CSR format  
    if hasattr(combined_adata.X, 'toarray'):
        if not sp.isspmatrix_csr(combined_adata.X):
            combined_adata.X = combined_adata.X.tocsr()
    
    print(f"Combined dataset shape: {combined_adata.shape}")
    print(f"Combined dataset value range: [{combined_adata.X.min():.3f}, {combined_adata.X.max():.3f}]")
    print(f"Cell types in combined dataset:")
    celltype_counts = combined_adata.obs[celltype_column].value_counts()
    for celltype, count in celltype_counts.items():
        print(f"  {celltype}: {count} cells")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the combined reference
    print(f"\nSaving combined reference to: {output_path}")
    combined_adata.write_h5ad(output_path, compression='gzip')
    
    print("Reference creation completed successfully!")
    return combined_adata

def apply_lognorm_tp10k(adata):
    """
    Apply log normalization with target sum 10,000 (lognorm tp10k).
    
    This is the standard normalization used by Seurat and scimilarity:
    1. Scale each cell to 10,000 total counts
    2. Add pseudocount of 1 
    3. Take natural logarithm
    
    Parameters:
    -----------
    adata : AnnData
        Input AnnData object with raw counts
        
    Returns:
    --------
    AnnData
        AnnData object with normalized data in .X
    """
    print("  Applying lognorm tp10k normalization...")
    
    # Ensure input is CSR format
    if hasattr(adata.X, 'toarray'):
        if not sp.isspmatrix_csr(adata.X):
            adata.X = adata.X.tocsr()
    
    # Calculate total counts per cell
    if hasattr(adata.X, 'toarray'):
        # Sparse matrix
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
    else:
        # Dense matrix
        total_counts = np.sum(adata.X, axis=1)
    
    # Avoid division by zero
    total_counts[total_counts == 0] = 1
    
    # Normalize to 10,000 counts per cell
    if hasattr(adata.X, 'toarray'):
        # For sparse matrices - work with CSR format
        normalized_X = adata.X.copy()
        # Normalize each row (cell) to sum to 10,000
        normalized_X = normalized_X.multiply(1/total_counts[:, None]).multiply(10000)
        # Ensure it's CSR format
        if not sp.isspmatrix_csr(normalized_X):
            normalized_X = normalized_X.tocsr()
    else:
        # For dense matrices
        normalized_X = (adata.X / total_counts[:, None]) * 10000
    
    # Add pseudocount and log transform: log(count + 1)
    if hasattr(normalized_X, 'toarray'):
        # For sparse matrices, we need to handle this carefully
        normalized_X.data = np.log1p(normalized_X.data)  # log1p is log(x + 1)
        # Ensure still CSR format after modification
        if not sp.isspmatrix_csr(normalized_X):
            normalized_X = normalized_X.tocsr()
    else:
        # For dense matrices
        normalized_X = np.log1p(normalized_X)
    
    # Update the AnnData object with CSR format
    adata.X = normalized_X
    
    print(f"  Normalization complete. Data range: [{adata.X.min():.3f}, {adata.X.max():.3f}]")
    print(f"  Matrix format: {type(adata.X)}")
    
    return adata

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a combined reference dataset from filtered adata files with celltypes')
    parser.add_argument('--celltype_column', type=str, default='celltype', 
                       help='Name of the obs column containing cell type information (default: celltype)')
    
    args = parser.parse_args()
    print(f"Using celltype column: {args.celltype_column}")
    
    create_reference_from_subset(celltype_column=args.celltype_column)
