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

def merge_counts(features_adata, counts_adata):
    #create two new columns in counts_adata.obs and set them to empty strings or nan
    counts_adata.obs['feature-10'] = ''
    counts_adata.obs['feature-4'] = ''
    #write the features_adata obs to counts_adata obs
    counts_adata.obs['feature-10'] = features_adata.obs['feature-10']
    counts_adata.obs['feature-4'] = features_adata.obs['feature-4']
    return counts_adata

def merge_features(features_adata, counts_adata):
    #subset the counts_adata to only include the obs that are in features_adata
    counts_adata = counts_adata[features_adata.obs_names]
    #copy all the obs columns from counts_adata to features_adata
    features_adata.obs = counts_adata.obs
    return features_adata

#print out all the env variables
for key, value in os.environ.items():
    print(f"{key}: {value}")

#find env features_dir and expressions_dir
features_dir_str = os.getenv('features_dirs')
expressions_dir_str = os.getenv('expression_dirs')

features_dir = ast.literal_eval(features_dir_str) if features_dir_str else []
expressions_dir = ast.literal_eval(expressions_dir_str) if expressions_dir_str else []

#print out the values and types
print(f"features_dir: {features_dir}, type: {type(features_dir)}")
print(f"expressions_dir: {expressions_dir}, type: {type(expressions_dir)}")
exit(0)

#find all files in features_dir
features_files = [f for f in os.listdir(features_dir) if os.path.isfile(os.path.join(features_dir, f))]


parser = argparse.ArgumentParser(description='Need to know where the base directory is for larry counts and the base directory for sc counts')

parser.add_argument('--upload-dir', required=True, type=str,
                    help='Specify the expression directory.')                    
parser.add_argument('--sc-dir', required=True, type=str,
                    help='Specify the sc directory.')
parser.add_argument('--overwrite', action='store_true', help='Overwrite existing files')    

# Parse the arguments
args = parser.parse_args()
                    
overwrite = False                  
sc_base_directory = args.sc_dir
upload_base_dir = args.upload_dir
overwrite = args.overwrite
methods_directories = ['kallisto','piscem/splicei','star']
methods_features_directories = ['kallisto','piscem.splicei','star']
#find all directries in sc_directory
sub_directories = [d for d in os.listdir(sc_base_directory) if os.path.isdir(os.path.join(sc_base_directory, d))]
for dir in sub_directories:
    for i in range(len(methods_directories)):
        method_features = methods_features_directories[i]
        counts_file = os.path.join(upload_base_dir, dir, method_features, 'unfiltered_counts.h5ad')
        print(f"Reading sc counts from {counts_file}")
        feature_file = os.path.join(upload_base_dir, dir, method_features, 'unfiltered_features.h5ad')
        print(f"Reading feature counts from {feature_file}")
        #make sure the files exist
        if not os.path.exists(counts_file):
            print(f"File {counts_file} does not exist")
            sys.exit(1)
        if not os.path.exists(feature_file):
            print(f"File {feature_file} does not exist")
            sys.exit(1)
        counts_adata = ad.read_h5ad(counts_file)
        feature_adata = ad.read_h5ad(feature_file)
        print(counts_adata)
        print(feature_adata)
        merged_features_adata=merge_features(feature_adata, counts_adata)
        print(merged_features_adata)
        #subset the adata object to only include rows where 'singlet_filtered' is True
        output_file= os.path.join(upload_base_dir, dir, method_features, 'unfiltered_features.h5ad')
        print(f"Writing merged features to {output_file}")
        merged_features_adata.write(output_file)
