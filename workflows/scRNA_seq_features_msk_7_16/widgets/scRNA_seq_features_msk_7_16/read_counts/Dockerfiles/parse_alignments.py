#!/usr/local/bin/python3
import os
import ast
import argparse
import scanpy as sc
import sys
import pyroe
import anndata as ad
import pandas as pd

def identify_experiment_type(dir):
    # check if the dir/counts_unfiltered exists and is a directory
    if os.path.exists(os.path.join(dir, 'counts_unfiltered')):
        if os.path.isdir(os.path.join(dir, 'counts_unfiltered')):
            return "kallisto"
    # check if dir/results/af_quant exists and is a directory
    if os.path.exists(os.path.join(dir, 'results/af_quant')):
        if os.path.isdir(os.path.join(dir, 'results/af_quant')):
            return "fry"
    #check if dir/Solo.out exists and is a directory
    if os.path.exists(os.path.join(dir, 'Solo.out')):
        if os.path.isdir(os.path.join(dir, 'Solo.out')):
            return "solo"
    #check if dir/outs exists and is a directory
    if os.path.exists(os.path.join(dir, 'outs')):
        if os.path.isdir(os.path.join(dir, 'outs')):
            return "cellranger"
    return None
    
def read_cellranger(directory):
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
    adata_raw.obs['filter'] = is_cell
    return adata_raw


def read_kallisto(directory):
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

def read_fry(frydir):    
    # Load the fry directory
    counts_directory = os.path.join(frydir, 'results/af_quant')
    custom_format = {'X' : ['U', 'S' ,'A'],
                      'spliced' : ['S'],
                      'unspliced' : ['U'],
                      'ambiguous' : ['A']
                    }
    fry_adata = pyroe.load_fry(counts_directory, output_format=custom_format)
    # Write the fry directory to the output path
    return fry_adata

def read_solo(directory):
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
    outputFeatureFilePath=os.path.join(inputFileDir, 'features.h5ad')
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
def process_directory(fulldir):
        fullpath = os.path.join(fulldir, 'counts.h5ad')  
        if os.path.exists(fullpath):
            if os.getenv("overwrite"):
                print(f"Overwriting the existing file {fullpath}")
            else:
                print(f"Skipping reading counts from  {fulldir} as the output file {fullpath} already exists.")
                return
        print(f"Processing directory {fulldir}")
        experiment_type = identify_experiment_type(fulldir)
        if experiment_type is None:
            print(f"Error: Could not identify the experiment type for the directory {fulldir}. Skipping this directory.")
            return
        print(f"Identified experiment type: {experiment_type}")
        print (f"Reading counts from {fulldir}")
        adata = readCounts(experiment_type, fulldir)
        print(f"counts data structure {adata}")
        #print first 5 var_names
        print(f"First 5 var_names {adata.var_names[:5]}")
        #print first 5 obs_names
        print(f"First 5 obs_names {adata.obs_names[:5]}")
        adata.write_h5ad(fullpath)
        print(f"Writing the counts to {fullpath}")
def readCounts(type, directory):
    readers = {
        'cellranger': read_cellranger,
        'kallisto': read_kallisto,
        'fry': read_fry,
        'solo': read_solo
    }
    
    if type not in readers:
        raise ValueError(f"Unsupported type: {type}")

    return readers[type](directory)    


alignsDir = os.getenv('alignsDir')

for dir in ast.literal_eval(alignsDir):
    #find the subdirectories in the directory
    subdirs=os.listdir(dir)
    for d in subdirs:
        #check if splicei or spliceu exits in the subdirectory
        if os.path.exists(os.path.join(dir, d, 'splicei')):
            print(f"Processing directory {d}")
            process_directory(os.path.join(dir, d, 'splicei'))
        elif os.path.exists(os.path.join(dir, d, 'spliceu')):
            print(f"Processing directory {d}")
            process_directory(os.path.join(dir, d, 'spliceu')) 
        else:
            process_directory(os.path.join(dir, d)) 
        
            


print("Processing complete.")

        
    


