#!/usr/local/bin/Rscript

library(reticulate)
library(Matrix)
library(DropletUtils)
library(scDblFinder)
library(Seurat)
library(anndata)
# Specify the path to the desired Python executable
use_python("/usr/bin/python3", required = TRUE)

# read the command line arguments which should be path of the h5ad file
args <- commandArgs(trailingOnly = TRUE)
h5ad_file <- args[1]
ad <- read_h5ad(h5ad_file)

counts <- as(t(ad$X), "dgCMatrix")
filtered_counts <- counts[, colSums(counts) > 0]

emptyDrops_result <- emptyDropsCellRanger(filtered_counts)
valid_cells <- !is.na(emptyDrops_result$FDR) & emptyDrops_result$FDR <= 0.01
non_empty_droplets <- filtered_counts[, valid_cells]
print( "The dimensions of the non-empty droplets are:" )
print( dim(non_empty_droplets) )
seurat_obj <- CreateSeuratObject(counts = non_empty_droplets)


#seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
sce <- as.SingleCellExperiment(seurat_obj)
sce <- scDblFinder(sce)
doublet_info <- sce$scDblFinder.class
counts_matrix <- seurat_obj@assays$RNA$counts
final_counts <- counts_matrix[, doublet_info == "singlet"]
final_counts <- as(final_counts, "dgCMatrix")
print( "The dimensions of the final counts matrix are:" )
print( dim(final_counts) )
# Write the valid cells to a text file
valid_barcodes <- colnames(final_counts)
#output path should be the same as the input path directory
output_path <- paste0(dirname(h5ad_file), "/filtered_barcodes.txt")
writeLines(valid_barcodes, output_path)