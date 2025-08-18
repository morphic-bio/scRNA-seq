import argparse
import numpy as np
import pandas as pd
import anndata as ad
import os
from scipy.sparse import issparse


def compare_anndata(adata1, adata2):
    differences = []

    # Compare shapes
    if adata1.shape != adata2.shape:
        differences.append(f"Shape mismatch: {adata1.shape} vs {adata2.shape}")

    # Compare X
    if not np.array_equal(adata1.X, adata2.X):
        differences.append("Difference in X matrix")

    # Compare obs
    if not adata1.obs.equals(adata2.obs):
        differences.append("Difference in obs DataFrame")

    # Compare var
    if not adata1.var.equals(adata2.var):
        differences.append("Difference in var DataFrame")

    # Compare uns
    if adata1.uns.keys() != adata2.uns.keys():
        differences.append("Difference in uns keys")
    else:
        for key in adata1.uns.keys():
            if not np.array_equal(adata1.uns[key], adata2.uns[key]):
                differences.append(f"Difference in uns[{key}]")

    # Compare obsm
    if adata1.obsm.keys() != adata2.obsm.keys():
        differences.append("Difference in obsm keys")
    else:
        for key in adata1.obsm.keys():
            if not np.array_equal(adata1.obsm[key], adata2.obsm[key]):
                differences.append(f"Difference in obsm[{key}]")

    # Compare varm
    if adata1.varm.keys() != adata2.varm.keys():
        differences.append("Difference in varm keys")
    else:
        for key in adata1.varm.keys():
            if not np.array_equal(adata1.varm[key], adata2.varm[key]):
                differences.append(f"Difference in varm[{key}]")

    # Compare layers
    if adata1.layers.keys() != adata2.layers.keys():
        differences.append("Difference in layers keys")
    else:
        for key in adata1.layers.keys():
            if not np.array_equal(adata1.layers[key], adata2.layers[key]):
                differences.append(f"Difference in layers[{key}]")

    # Report differences
    if not differences:
        print("The AnnData objects are identical.")
    else:
        print("Differences found between the AnnData objects:")
        for difference in differences:
            print(difference)




def ensure_string_dtype(obs):
    """
    Ensures all columns in the obs DataFrame are of appropriate string types if needed.
    """
    for col in obs.columns:
        if obs[col].dtype == 'object' and not pd.api.types.is_string_dtype(obs[col]):
            obs[col] = obs[col].astype(str)
    return obs

# Parse arguments to get input file
parser = argparse.ArgumentParser(description='Create a filtered AnnData object')
parser.add_argument('--input_file', type=str, help='Input AnnData file', required=True)
args = parser.parse_args()

input_file = args.input_file
directory = os.path.dirname(input_file)

# Read the AnnData object from input file
#print(f"Reading the AnnData object from {input_file}")
adata = ad.read_h5ad(input_file)
print(adata)

# Read the features AnnData object
features_file = os.path.join(directory, 'merged_features.h5ad')
print(f"Reading the features AnnData object from {features_file}")
features_adata = ad.read_h5ad(features_file)
print(features_adata)
exit(0)
#extract the layers from the adata object
layer_data = features_adata.layers['sc_layer']
layer_adata = ad.AnnData(X=layer_data, obs=features_adata.obs, var=features_adata.var)
print(layer_adata)
del features_adata.layers['sc_layer']
raw_features_adata = ad.read_h5ad(os.path.join(directory, 'raw_features.h5ad'))
#print(layer_adata.obs)
#print(layer_adata.var)
compare_anndata(features_adata, raw_features_adata)
# Merge the obs from adata into features_adata
for key in adata.obs.keys():
    raw_features_adata.obs[key] = adata.obs[key]
    features_adata.obs[key] = adata.obs[key]


# Save the merged AnnData object
output_file = os.path.join(directory, 'unfiltered_features.h5ad')
raw_features_adata.write(output_file)
print(f"Wrote raw AnnData object to {output_file}")
#features_adata.write(output_file)
#print(f"Wrote features AnnData object to {output_file}")
exit(0)
# Create a binary filter from adata.obs['singlet_filtered'] and use it to subset features_adata
singlet_filtered = adata.obs['singlet_filtered']
features_adata_filtered = features_adata[singlet_filtered].copy()

# Ensure all obs columns are of appropriate types for the filtered data
features_adata_filtered.obs = ensure_string_dtype(features_adata_filtered.obs)
print(f"Layers in adata: {adata.layers.keys()}")
print(f"Layers in features_adata: {features_adata.layers.keys()}")
# Save the filtered AnnData object
output_file = os.path.join(directory, 'filtered_features.h5ad')
print(f"Writing filtered AnnData object to {output_file}")
features_adata_filtered.write(output_file)
