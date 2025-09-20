#!/usr/bin/env python3
import os
import ast
import argparse
import scanpy as sc
import sys
import pyroe
import anndata as ad
import pandas as pd
import multiprocessing
from itertools import product

def identify_experiment_type(dir):
    if os.path.exists(os.path.join(dir) and os.path.isdir(os.path.basename(dir)) and os.path.isdir(os.path.join(dir, 'demux'))):
        return "flex"
    # check if the dir/counts_unfiltered exists and is a directory
    if os.path.exists(os.path.join(dir, 'counts_unfiltered')):
        if os.path.isdir(os.path.join(dir, 'counts_unfiltered')):
            return "kallisto"
    # check if dir/results/af_quant exists and is a directory
    if os.path.exists(os.path.join(dir, 'af_quant')):
        if os.path.isdir(os.path.join(dir, 'af_quant')):
            return "fry"
    #check if dir/Solo.out exists and is a directory
    if os.path.exists(os.path.join(dir, 'Solo.out')):
        if os.path.isdir(os.path.join(dir, 'Solo.out')):
            return "solo"
    #check if dir/outs exists and is a directory
    if os.path.exists(os.path.join(dir, 'outs')):
        if os.path.isdir(os.path.join(dir, 'outs')):
            return "cellranger"
    print(f"Error: Could not identify the experiment type for the directory {dir}. Skipping this directory.")
    return None
def read_cellranger(directory):
    print("Reading cellranger data")
    filtered_dir = os.path.join(directory, 'outs/filtered_feature_bc_matrix')
    raw_dir = os.path.join(directory, 'outs/raw_feature_bc_matrix')
    # check if these directories exist - if not return an error message
    if not os.path.exists(filtered_dir):
        print(f"Error: The directory '{filtered_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(raw_dir): 
        print(f"Error: The directory '{raw_dir}' does not exist.", file=sys.stderr)
        sys.exit(1) 
    
    # Read filtered data
    adata_filtered = sc.read_10x_mtx(filtered_dir, var_names='gene_ids', cache=True)
    # Read raw data
    adata_raw = sc.read_10x_mtx(raw_dir, var_names='gene_ids', cache=True)
    # Create a mask that indicates whether the barcode is a cell or not
    cell_barcodes = adata_filtered.obs_names
    is_cell = adata_raw.obs_names.isin(cell_barcodes)
    # Add the mask to the raw data AnnData object
    adata_raw.obs['is_cell'] = is_cell
    return adata_raw


def read_kallisto(directory):
    print ("Reading kallisto data")
    counts_directory = os.path.join(directory, 'counts_unfiltered') 
    total_file=os.path.join(counts_directory, 'cells_x_genes.total.mtx')
    unspliced_file=os.path.join(counts_directory, 'cells_x_genes.nascent.mtx')
    spliced_file=os.path.join(counts_directory, 'cells_x_genes.mature.mtx')
    kallisto_ambiguous_file=os.path.join(counts_directory, 'cells_x_genes.ambiguous.mtx')
    bpath=os.path.join(counts_directory, 'cells_x_genes.barcodes.txt')
    gpath=os.path.join(counts_directory, 'cells_x_genes.genes.txt')
    #read the matrices in
    x = sc.read_mtx(total_file)
    u = sc.read_mtx(unspliced_file)
    s = sc.read_mtx(spliced_file)
    a = sc.read_mtx(kallisto_ambiguous_file)
    #read the gene and barcode files in
    bc = pd.read_csv(bpath, names=['barcodes'])
    g = pd.read_csv(gpath, names=['gene_ids'])
    g.gene_ids = g.gene_ids.str.split('.').str.get(0)
    gene_names=g.gene_ids.to_list()
    barcode_names=bc.barcodes.to_list()
    
    
    bc = bc.set_index('barcodes')
    g = g.set_index('gene_ids')
    
    x.var = g
    x.obs = bc
    adata = ad.AnnData(X=x.X)
    adata.var_names = gene_names
    adata.obs_names = barcode_names
    adata.layers['unspliced'] = u.X
    adata.layers['spliced'] = s.X
    adata.layers['ambiguous'] = a.X
    return adata
    print("Reading salmon data")
    alevin_dir = os.path.join(directory, 'alevin')
    mtx_file = os.path.join(alevin_dir, 'quants_mat.mtx')
    genes_file = os.path.join(alevin_dir, 'quants_mat_rows.txt')
    barcodes_file = os.path.join(alevin_dir, 'quants_mat_cols.txt')

    adata_matrix = sc.read_mtx(mtx_file)  # genes x cells
    genes = pd.read_csv(genes_file, header=None, sep='\t', names=['gene_ids'])
    barcodes = pd.read_csv(barcodes_file, header=None, sep='\t', names=['barcodes'])

    adata = ad.AnnData(adata_matrix.T)  # cells x genes
    adata.var_names = genes['gene_ids'].values
    adata.obs_names = barcodes['barcodes'].values

    # Salmon alevin doesn't output spliced/unspliced layers in this default format.
    return adata

def read_fry(frydir):
    print("Reading fry data")
    # Load the fry directory
    counts_directory = os.path.join(frydir, 'af_quant')
    custom_format = {'X' : ['U', 'S' ,'A'],
                      'spliced' : ['S'],
                      'unspliced' : ['U'],
                      'ambiguous' : ['A']
                    }
    fry_adata = pyroe.load_fry(counts_directory, output_format=custom_format)
    # Write the fry directory to the output path
    return fry_adata

def read_solo(directory):
    print("Reading solo data")
    counts_directory = os.path.join(directory, 'Solo.out')
    total_raw_file=os.path.join(counts_directory, 'GeneFull/raw/matrix.mtx')
    barcodes_file=os.path.join(counts_directory, 'GeneFull/raw/barcodes.tsv')
    genes_file=os.path.join(counts_directory, 'GeneFull/raw/features.tsv')
    
    adata_matrix= sc.read_mtx(total_raw_file)
    # Optionally, read the gene names and cell barcodes if available
    genes = pd.read_csv(genes_file, header=None, sep='\t')
    barcodes = pd.read_csv(barcodes_file, header=None, sep='\t')

    # Create an AnnData object
    adata = sc.AnnData(adata_matrix.T)  # Transpose to make genes rows and cells columns
    adata.var_names = genes[0]
    adata.obs_names = barcodes[0]
    
    #create mask to filter out the cells that are not in the filtered matrix
    cell_barcodes = pd.read_csv(os.path.join(counts_directory, 'GeneFull/filtered/barcodes.tsv'), header=None, sep='\t')[0]
    is_cell = adata.obs_names.isin(cell_barcodes)
    adata.obs['filter'] = is_cell
    
    #add layers to the adata object
    spliced_file=os.path.join(counts_directory, 'Velocyto/raw/spliced.mtx')
    unspliced_file=os.path.join(counts_directory, 'Velocyto/raw/unspliced.mtx')
    ambiguous_file=os.path.join(counts_directory, 'Velocyto/raw/ambiguous.mtx')
    
    spliced = sc.read_mtx(spliced_file)
    unspliced = sc.read_mtx(unspliced_file)
    ambiguous = sc.read_mtx(ambiguous_file)
    
    #add layers to the adata object
    adata.layers['spliced'] = spliced.X.T
    adata.layers['unspliced'] = unspliced.X.T
    adata.layers['ambiguous'] = ambiguous.X.T
    return adata

def filter_features(adata,inputFilePath):
    inputFileDir=os.path.dirname(inputFilePath)
    outputFeatureFilePath=os.path.join(inputFileDir, 'features_gex.h5ad')
    print("filtering the AnnData object")
    #find the rows do not contain 'ENSG' and save them to adata as adata.var['features'] 
    feature_names = ~adata.var_names.str.contains('ENSG')
    #check if this is empty
    if feature_names.sum() == 0:
        print("No rows contain 'ENSG' in the label", file=sys.stderr)
        return adata
    #subset the against feature_names and save to an new AnnData object
    features_adata = adata[:, feature_names]
    #remove the layers from the AnnData object
    features_adata.layers = {}
    print(features_adata)
    #write the AnnData object to a h5ad file
    features_adata.write_h5ad(outputFeatureFilePath)
    #subset the AnnData object  to only include the rows that contain 'ENSG' in the label and subset the layers     adata.layers['unsplice'], adata.layers['spliced'], adata.layers['ambiguous'] to only include the rows that contain 'ENSG' in the label 
    return adata[:, ~feature_names]

def process_discovered_directory(args):
    base_dir, sub_dir = args
    path_splicei = os.path.join(base_dir, sub_dir, 'splicei')
    path_spliceu = os.path.join(base_dir, sub_dir, 'spliceu')
    found_splice_dir = False

    if os.path.exists(path_splicei):
        print(f"Processing directory {sub_dir}")
        process_directory(path_splicei)
        found_splice_dir = True
    
    if os.path.exists(path_spliceu):
        print(f"Processing directory {sub_dir}")
        process_directory(path_spliceu)
        found_splice_dir = True

    if not found_splice_dir:
        process_directory(os.path.join(base_dir, sub_dir))

def process_directory(fulldir):
        fullpath = os.path.join(fulldir, 'unfiltered_counts.h5ad') 
        if os.path.exists(fullpath):
            if os.getenv("overwrite", "false").lower() in ('true', '1', 't'):
                print(f"Overwriting the existing file {fullpath}")
            else:
                print(f"Skipping reading counts from  {fulldir} as the output file {fullpath} already exists.")
                return
        print(f"Processing directory {fulldir}",flush=True)
        experiment_type = identify_experiment_type(fulldir)
        print(f"Experiment type: {experiment_type}")
        if experiment_type is None:
            print(f"Error: Could not identify the experiment type for the directory {fulldir}. Skipping this directory.")
            return
        print(f"Identified experiment type: {experiment_type}")
        if (experiment_type == 'flex'):
            gex_adata, demux_adata = read_flex(fulldir)
            print(f"Writing flex gex data to {fullpath}")
            gex_adata.write_h5ad(os.path.join(fulldir, 'gex.h5ad'))
            print(f"Writing flex demux data to {fullpath}")
            demux_adata.write_h5ad(os.path.join(fulldir, 'demux.h5ad'))
            return
        else:
            adata = readCounts(experiment_type, fulldir)

        filter_features_flag = os.getenv('filter_features', 'false').lower() in ('true', '1', 't')
        if filter_features_flag:
            print(f"Filtering features for {fulldir}")
            adata = filter_features(adata, fullpath)

        print(f"Writing counts to {fullpath}")
        adata.write_h5ad(fullpath)
def read_gex_subdirectory(directory):
    print("Reading flex gex data")
    print(f"Reading flex gex data from {directory}")
    adata = sc.read_mtx(os.path.join(directory, 'matrix.mtx'))
    adata = adata.T
    print(f"Reading genes from {os.path.join(directory, 'features.tsv')}")
    genes = pd.read_csv(os.path.join(directory, 'features.tsv'), header=None, sep='\t', names=['gene_ids'])
    barcodes = pd.read_csv(os.path.join(directory, 'barcodes.tsv'), header=None, sep='\t', names=['barcodes'])
    adata.var_names = genes['gene_ids'].values
    adata.obs_names = barcodes['barcodes'].values
    return adata

def read_demux_subdirectory(directory):
    print("Reading flex demux data")
    print(f"Reading flex demux data from {directory}")
    adata = sc.read_mtx(os.path.join(directory, 'matrix.mtx'))
    print(f"Reading sample_ids from {os.path.join(directory, 'features.tsv')}")
    sample_ids = pd.read_csv(os.path.join(directory, 'features.tsv'), header=None, sep='\t', names=['sample_ids'])
    print(f"Reading barcodes from {os.path.join(directory, 'barcodes.tsv')}")
    barcodes = pd.read_csv(os.path.join(directory, 'barcodes.tsv'), header=None, sep='\t', names=['barcodes'])
    
    print(f"Matrix shape: {adata.shape}")
    print(f"Number of sample_ids: {len(sample_ids)}")
    print(f"Number of barcodes: {len(barcodes)}")
    adata = adata.T
    
    adata.var_names = sample_ids['sample_ids'].values
    adata.obs_names = barcodes['barcodes'].values
    return adata
def read_flex(directory):
    demux_directory = os.path.join(directory, 'demux')
    gex_directory = os.path.join(directory, 'gex')
    print(f"Reading flex demux data from {demux_directory}")
    demux_adata = read_demux_subdirectory(demux_directory)
    print(f"Reading flex gex data from {gex_directory}")
    gex_adata = read_gex_subdirectory(gex_directory)
    return gex_adata, demux_adata

def readCounts(type, directory):
    readers = {
        'cellranger': read_cellranger,
        'kallisto': read_kallisto,
        'fry': read_fry,
        'solo': read_solo,
    }
    if type not in readers:
        raise ValueError(f"Unsupported type: {type}")
    print(f"Reading counts from {directory}",flush=True) 
    return readers[type](directory)    


alignsdir = os.getenv('alignsdir')
if '[' in alignsdir:
    dirs = ast.literal_eval(alignsdir)
else: 
    dirs = [alignsdir]
print(f"Processing directories {dirs}")

subdirs_to_process = []
for dir_path in dirs:
    print(f"Scanning directory {dir_path}")
    # find the subdirectories in the directory
    subdirs = os.listdir(dir_path)
    subdirs_to_process.extend(product([dir_path], subdirs))

try:
    # Get nThreads from environment variable, default to number of CPUs
    nThreads = int(os.getenv('nThreads', multiprocessing.cpu_count()))
except (ValueError, TypeError):
    nThreads = multiprocessing.cpu_count()

# Ensure nThreads is not greater than the number of tasks
if nThreads > len(subdirs_to_process):
    nThreads = len(subdirs_to_process)

if nThreads > 0:
    print(f"Processing {len(subdirs_to_process)} subdirectories in parallel with {nThreads} threads.")
    with multiprocessing.Pool(processes=nThreads) as pool:
        pool.map(process_discovered_directory, subdirs_to_process)
else:
    print("No subdirectories to process.")

print("Processing complete.")



