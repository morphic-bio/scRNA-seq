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
output_path <- paste0(dirname(h5ad_file), "/non_empty_barcodes.txt")
writeLines(colnames(non_empty_droplets), output_path)

seurat_obj <- CreateSeuratObject(counts = non_empty_droplets)
#seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
sce <- as.SingleCellExperiment(seurat_obj)
sce <- scDblFinder(sce)
doublet_info <- sce$scDblFinder.class
doublet_scores <- sce$scDblFinder.score

# Create a data frame to store the barcode, classification, and score
barcodes <- colnames(sce)
doublet_results <- data.frame(
  Barcode = barcodes,
  Classification = doublet_info,
  Score = doublet_scores
)
# Write the doublet information (classification and score) to a file
output_path_with_scores <- paste0(dirname(h5ad_file), "/filtered_barcodes_with_scores.txt")
write.table(doublet_results, output_path_with_scores, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

counts_matrix <- seurat_obj@assays$RNA$counts
#find doublets
doublet_barcodes <- colnames(counts_matrix)[doublet_info == "doublet"]
# Write the doublet barcodes to a text file
output_path <- paste0(dirname(h5ad_file), "/doublet_barcodes.txt")
writeLines(doublet_barcodes, output_path)