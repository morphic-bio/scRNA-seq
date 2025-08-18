import pandas as pd
import os
import anndata as ad
import scipy.sparse as sp

def load_msk_village_celltypes(csv_path='/mnt/pikachu/scRNA-seq/data/MSK_KO_village_meta_data.csv'):
    """
    Load MSK village cell type data from CSV file.
    
    Parameters:
    -----------
    csv_path : str
        Path to the CSV file containing MSK village metadata
        
    Returns:
    --------
    dict
        Nested dictionary where celltypes[sample_name][barcode] = celltype
    """
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    # Initialize the nested dictionary
    celltypes = {}
    
    # Process each row
    for _, row in df.iterrows():
        # Extract and clean the fields (remove quotes)
        sample_name = str(row['sample_name']).strip('"')
        full_barcode = str(row['barcode']).strip('"')
        celltype = str(row['celltype']).strip('"')
        
        # Extract the ACGT string from barcode (between last underscore and hyphen)
        # Example: "Sample_B_WT_AAACCCAAGGAAGTAG-1" -> "AAACCCAAGGAAGTAG"
        if '_' in full_barcode and '-' in full_barcode:
            # Find the last underscore
            last_underscore_idx = full_barcode.rfind('_')
            # Find the hyphen after the last underscore
            hyphen_idx = full_barcode.find('-', last_underscore_idx)
            
            if last_underscore_idx != -1 and hyphen_idx != -1:
                barcode = full_barcode[last_underscore_idx + 1:hyphen_idx]
            else:
                # Fallback if format is unexpected
                barcode = full_barcode
        else:
            # Fallback if format is unexpected
            barcode = full_barcode
        
        # Initialize sample_name dictionary if it doesn't exist
        if sample_name not in celltypes:
            celltypes[sample_name] = {}
        
        # Store the celltype
        celltypes[sample_name][barcode] = celltype
    
    print(f"Loaded cell types for {len(celltypes)} samples")
    for sample, barcodes in celltypes.items():
        print(f"  {sample}: {len(barcodes)} cells")
    
    return celltypes

def ensure_csr_format(adata):
    """Ensure AnnData X matrix and all layers are in CSR format."""
    # Fix main matrix
    if hasattr(adata.X, 'toarray') and not sp.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()
    
    # Fix layers if they exist
    if hasattr(adata, 'layers'):
        for layer_name in adata.layers:
            layer_matrix = adata.layers[layer_name]
            if hasattr(layer_matrix, 'toarray') and not sp.isspmatrix_csr(layer_matrix):
                adata.layers[layer_name] = layer_matrix.tocsr()
    
    return adata

# Example usage
if __name__ == "__main__":
    celltypes = load_msk_village_celltypes()
    #samples are the keys of the celltypes dictionary
    samples = list(celltypes.keys())
    alignments_dir = os.getenv('alignments_dir')
    if not alignments_dir:
        alignments_dir = '/mnt/pikachu/storage/scRNAseq_output/Alignments'
    for sample in samples:
        alignments_sample_dir = os.path.join(alignments_dir, sample)
        methods = [d for d in os.listdir(alignments_sample_dir) if os.path.isdir(os.path.join(alignments_sample_dir, d))]
        for method in methods:
            counts_h5ad = os.path.join(alignments_sample_dir, method, 'merged_counts.h5ad')
            counts_adata = ad.read_h5ad(counts_h5ad)
            
            # Find intersection of barcodes between dictionary and adata
            adata_barcodes = counts_adata.obs_names
            dict_barcodes = set(celltypes[sample].keys())
            
            # Create celltype column initialized with empty string or NaN
            counts_adata.obs['user_celltype'] = ''
            
            # Map celltypes for matching barcodes
            for barcode in adata_barcodes:
                if barcode in dict_barcodes:
                    counts_adata.obs.loc[barcode, 'user_celltype'] = celltypes[sample][barcode]
            
            # Print matching statistics
            matching_barcodes = len(set(adata_barcodes).intersection(dict_barcodes))
            print(f"Sample {sample}, method {method}: {matching_barcodes}/{len(adata_barcodes)} barcodes matched")
            
            # Ensure CSR format before writing to prevent AnnData warnings
            counts_adata = ensure_csr_format(counts_adata)
            counts_adata.write(os.path.join(alignments_sample_dir, method, 'merged_counts_with_celltypes.h5ad'))