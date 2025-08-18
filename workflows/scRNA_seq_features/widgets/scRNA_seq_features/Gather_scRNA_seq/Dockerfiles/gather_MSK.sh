#!/bin/bash

MSK_dir=/mnt/pikachu/MSK_village
samples=( A B C D E F G_1 G_2 H I J L_1 L_2 )
alignments_dir=/mnt/pikachu/storage/scRNAseq_output/Alignments
features_dir=/mnt/pikachu/processing/morphic-test/test/counts
ca_assessment_dir=/mnt/pikachu/storage/scRNAseq_output/batch_analysis
mkdir -p $MSK_dir/QC/cell_assignments
echo "cp $ca_assessment_dir/summary* $MSK_dir/QC/cell_assignments"

for sample in ${samples[@]}; do
    mkdir -p $MSK_dir/$sample/QC/features
    mkdir -p $MSK_dir/$sample/QC/cell_assignments
    mkdir -p $MSK_dir/$sample/QC/counts
    echo "cp $alignments_dir/$sample/star/merged_counts_with_celltypes_classified.h5ad $MSK_dir/$sample/counts.h5ad"
    cp $alignments_dir/$sample/star/merged_counts_with_celltypes_classified.h5ad $MSK_dir/$sample/counts.h5ad
    echo "cp $alignments_dir/$sample/star/merged_features.h5ad $MSK_dir/$sample/features.h5ad"
    cp $alignments_dir/$sample/star/merged_features.h5ad $MSK_dir/$sample/features.h5ad
    echo "cp $alignments_dir/$sample/star/gene_quantile* $MSK_dir/$sample/QC/counts"
    cp $alignments_dir/$sample/star/gene_quantile* $MSK_dir/$sample/QC/counts
    echo "cp $features_dir/$sample/heatmap.png $MSK_dir/$sample/QC/features"
    cp $features_dir/$sample/heatmap.png $MSK_dir/$sample/QC/features
    echo "cp $features_dir/$sample/feature_sequences.txt $MSK_dir/$sample/QC/features"
    cp $features_dir/$sample/stats.txt $MSK_dir/$sample/QC/features
    echo "cp $features_dir/$sample/feature_sequences.txt $MSK_dir/$sample/QC/features"
    cp $features_dir/$sample/feature_sequences.txt $MSK_dir/$sample/QC/features
    echo "cp $features_dir/$sample/feature_sequences.txt $MSK_dir/$sample/QC/features"
    cp $features_dir/$sample/feature_sequences.txt $MSK_dir/$sample/QC/features
    echo "cp $ca_assessment_dir/${sample}_star/* $MSK_dir/$sample/QC/cell_assignments"
    cp $ca_assessment_dir/${sample}_star/* $MSK_dir/$sample/QC/cell_assignments
done
