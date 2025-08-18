import argparse
import numpy as np
import pandas as pd
import anndata as ad
import os
from scipy.sparse import issparse
from scipy.sparse import csr_matrix


def align_obs_with_zeros(adata, features_adata):
    """
    Aligns the observations (obs) of features_adata with adata. Adds zeros for missing values and removes extra obs.

    Parameters:
    adata (AnnData): The reference AnnData object.
    features_adata (AnnData): The AnnData object to be aligned.

    Returns:
    AnnData: The aligned AnnData object.
    """
    common_obs = adata.obs_names.intersection(features_adata.obs_names)
    missing_obs = adata.obs_names.difference(features_adata.obs_names)
    
    print (f"Common obs: {len(common_obs)}")
    print (f"Missing obs: {len(missing_obs)}")
    
    
    # Subset features_adata to only include common_obs
    features_adata = features_adata[common_obs, :]
    if len(missing_obs) == 0:
        return features_adata
    # Create a zero matrix for missing observations]
    if issparse(features_adata.X):
        zero_data = csr_matrix((len(missing_obs), features_adata.shape[1]), dtype=features_adata.X.dtype)
    else:
        zero_data = np.zeros((len(missing_obs), features_adata.shape[1]), dtype=features_adata.X.dtype)

    # Create a new AnnData object with zeros for missing observations
    missing_adata = ad.AnnData(X=zero_data, obs=pd.DataFrame(index=missing_obs), var=features_adata.var)

    # Concatenate the two AnnData objects
    aligned_adata = features_adata.concatenate(missing_adata)

    # Ensure the final obs are in the same order as in adata
    aligned_adata = aligned_adata[adata.obs_names, :]

    return aligned_adata


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

# Read the features AnnData object
features_file = os.path.join(directory, 'merged_features.h5ad')


print(f"Reading the features AnnData object from {features_file}")
features_adata = ad.read_h5ad(features_file)
print(features_adata)      
features_adata=align_obs_with_zeros(adata, features_adata)
print(features_adata)
for key in adata.obs.keys():
    features_adata.obs[key] = adata.obs[key]
# Save the merged AnnData object
output_file = os.path.join(directory, 'unfiltered_features.h5ad')
features_adata.write(output_file)
print(f"Wrote raw AnnData object to {output_file}")
# Create a binary filter from adata.obs['singlet_filtered'] and use it to subset features_adata]
singlet=adata.obs['singlet']
singlet_filtered = adata.obs['singlet_filtered']
features_adata_filtered = features_adata[singlet].copy()
# Save the filtered AnnData object
output_file = os.path.join(directory, 'filtered_features.h5ad')
print(f"Writing filtered AnnData object to {output_file}")
features_adata_filtered.write(output_file)
output_file = os.path.join(directory, 'filtered_hq_features.h5ad')
features_adata_filtered = features_adata[singlet_filtered].copy()
print(f"Writing filtered AnnData object to {output_file}")
features_adata_filtered.write(output_file)