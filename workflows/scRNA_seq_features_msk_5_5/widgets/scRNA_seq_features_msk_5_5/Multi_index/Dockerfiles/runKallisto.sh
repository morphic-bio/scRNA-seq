#!/bin/bash
export AF_SAMPLE_DIR=/data/af_test_workdir
export ALEVIN_FRY_HOME="$AF_SAMPLE_DIR/af_home"
mkdir -p $AF_SAMPLE_DIR
mkdir -p $ALEVIN_FRY_HOME
simpleaf set-paths

cd $AF_SAMPLE_DIR
FASTQ_DIR="/data/fastqs/pbmc_1k_v3_fastqs"


REF_DIR="/data/run_cellranger_count/refdata-gex-GRCh38-2020-A/"
IDX_DIR_SPLICEI="$AF_SAMPLE_DIR/human-2020-A_splicei"
IDX_DIR_SPLICEU="$AF_SAMPLE_DIR/human-2020-A_spliceu"
simpleaf index \
--output $IDX_DIR_SPLICEI \
--fasta $REF_DIR/fasta/genome.fa \
--gtf $REF_DIR/genes/genes.gtf \
--rlen 91 \
--threads 16 \
--use-piscem  # remove this if missing piscem

simpleaf index \
--output $IDX_DIR_SPLICEU \
--fasta $REF_DIR/fasta/genome.fa \
--gtf $REF_DIR/genes/genes.gtf \
--threads 16 \
--ref-type spliceu \
--use-piscem  # remove this if missing piscem
