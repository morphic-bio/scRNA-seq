#!/bin/bash

UCSF_dir=/mnt/pikachu/UCSF-perturb
samples=( $(ls /storage/features) )
alignments_dir=/storage/scRNAseq_output/Alignments
features_dir=/storage/features

for sample in ${samples[@]}; do
    mkdir -p $UCSF_dir/$sample/QC/features
    mkdir -p $UCSF_dir/$sample/QC/cell_assignments
    mkdir -p $UCSF_dir/$sample/QC/counts
    echo "cp $alignments_dir/$sample/star/merged_counts.h5ad $UCSF_dir/$sample/counts.h5ad"
    cp $alignments_dir/$sample/star/merged_counts.h5ad $UCSF_dir/$sample/counts.h5ad
    echo "cp $alignments_dir/$sample/star/merged_features.h5ad $UCSF_dir/$sample/features.h5ad"
    cp $alignments_dir/$sample/star/merged_features.h5ad $UCSF_dir/$sample/features.h5ad
    echo "cp $alignments_dir/$sample/star/gene_quantile* $UCSF_dir/$sample/QC/counts"
    cp $alignments_dir/$sample/star/gene_quantile* $UCSF_dir/$sample/QC/counts
    echo "cp $features_dir/$sample/stats.txt $UCSF_dir/$sample/QC/features"
    cp $features_dir/$sample/stats.txt $UCSF_dir/$sample/QC/features
    echo "cp $features_dir/$sample/heatmap.png $UCSF_dir/$sample/QC/features"
    cp $features_dir/$sample/heatmap.png $UCSF_dir/$sample/QC/features
    echo "cp $features_dir/$sample/feature_sequences.txt $UCSF_dir/$sample/QC/features"
    cp $features_dir/$sample/feature_sequences.txt $UCSF_dir/$sample/QC/features
done
