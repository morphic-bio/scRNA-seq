#!/usr/bin/env python3
import argparse
#note that scanpy seems to be masked somehow and does not show up so we use anndata instead
import anndata as ad
import os
import numpy as np
import cellbender.remove_background.downstream as cb
import scipy.sparse as sp


def add_denoised_layer(adata, cellbender_file, layer_name='denoised'):  
    # Load CellBender output as AnnData object
    cb_adata = cb.anndata_from_h5(cellbender_file)
    denoised_matrix = cb_adata.X  # This can be a sparse matrix

    # Print the type and shape for debugging
    print(f"type of denoised_matrix: {type(denoised_matrix)}")
    print(f"shape of denoised_matrix: {denoised_matrix.shape}")

    # Map observation names to indices
    obs_name_to_idx = {name: idx for idx, name in enumerate(adata.obs_names)}
    
    # Get the corresponding integer indices for the subset observations
    subset_indices = np.array([obs_name_to_idx[name] for name in cb_adata.obs_names])

    # Create a sparse matrix with the denoised values at the correct indices
    rows, cols = denoised_matrix.nonzero()
    data = denoised_matrix.data

    # Adjust the row indices based on the mapping to the original `adata`
    adjusted_rows = subset_indices[rows]

    # Create the new sparse matrix with values only in the correct positions
    new_layer = sp.csr_matrix((data, (adjusted_rows, cols)), shape=adata.shape)

    # Print the shape of the new layer
    print(f"new_layer shape: {new_layer.shape}")
    
    # Add the new layer to the AnnData object
    adata.layers[layer_name] = new_layer
    
    return adata



def main(raw_counts_h5ad, cellbender_file, output_h5ad, layer_name):
    # Check if the output file already exists
    if not os.getenv("overwrite_layer") and os.path.exists(output_h5ad):
        print(f"The output file {output_h5ad} already exists.")
        print (f"Checking if the layer {layer_name} already exists in the AnnData object.")
        adata = ad.read_h5ad(output_h5ad)
        if layer_name in adata.layers.keys():
            print(f"Skipping adding the layer {layer_name} as it already exists.")
            return
    # Load the old AnnData object
    adata = ad.read_h5ad(raw_counts_h5ad)
    # Add the denoised layer
    adata = add_denoised_layer(adata, cellbender_file, layer_name)
    
    # Save the updated AnnData object to a new file
    adata.write(output_h5ad)
    print(f"Updated AnnData object saved to {output_h5ad}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add a denoised layer to an AnnData object using CellBender results.")
    parser.add_argument('-i', required=True, help="Path to the raw counts AnnData (.h5ad) file.")
    parser.add_argument('-c', required=True, help="Path to the CellBender denoised data (.h5) file.")
    parser.add_argument('-o', required=True, help="Path to save the updated AnnData object (.h5ad) file.")
    parser.add_argument('-l', required=True, help="Name of the layer to add.")
    
    args = parser.parse_args()
    print(args)
    main(args.i, args.c, args.o, args.l)
