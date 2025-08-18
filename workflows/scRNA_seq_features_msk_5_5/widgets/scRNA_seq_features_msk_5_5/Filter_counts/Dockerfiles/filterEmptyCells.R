#!/usr/local/bin/Rscript
#this was written to fix the previous version of the script where both empty cell and doublet removal were done in one step
library(reticulate)
library(Matrix)
library(DropletUtils)
library(anndata)
# Specify the path to the desired Python executable
use_python("/usr/bin/python3", required = TRUE)

# read the command line arguments which should be path of the h5ad file
args <- commandArgs(trailingOnly = TRUE)
h5ad_file <- args[1]
ad <- read_h5ad(h5ad_file)
if ( class(ad$X)[1] != "dgCMatrix" ) {
  counts <- as(t(ad$X), "CsparseMatrix")
  counts <- as(counts, "dgCMatrix")
}else{
  counts <- as(t(ad$X), "dgCMatrix")
}
filtered_counts <- counts[, colSums(counts) > 0]

emptyDrops_result <- emptyDropsCellRanger(filtered_counts)
valid_cells <- !is.na(emptyDrops_result$FDR) & emptyDrops_result$FDR <= 0.01
non_empty_droplets <- filtered_counts[, valid_cells]
print( "The dimensions of the non-empty droplets are:" )
print( dim(non_empty_droplets) )
non_empty_barcodes <- colnames(non_empty_droplets)
output_path <- paste0(dirname(h5ad_file), "/non_empty_barcodes.txt")
writeLines(non_empty_barcodes, output_path)
#read in the valid barcodes from the previous step

#fixes previous version where both empty cell and doublet removal were done in one step
input_path <- paste0(dirname(h5ad_file), "/filtered_barcodes.txt")
if (file.exists(input_path)) {
  valid_barcodes <- readLines(input_path)
  #find subset of non_empty_barcodes that are not in valid_barcodes
  doublet_barcodes <- non_empty_barcodes[!non_empty_barcodes %in% valid_barcodes]
  output_path <- paste0(dirname(h5ad_file), "/doublet_barcodes.txt")
  writeLines(doublet_barcodes, output_path)
}
#exit the script

