#!/bin/bash

# raw inputs

#echo all the environment variables
echo "environment variables"
env

#echo all the arguments
echo "arguments"
echo "$@"


#find all the environment variables that are directories and make them if they are reasonable strings

#set the gencode version
gencode=44

#set the ensembl release
ens_release=114

[ -z "$work_dir" ] && work_dir="/data/scRNAseq_output"
indices_dir="$work_dir/indices"
download_dir="$work_dir/genome"
raw_fa="$download_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
raw_gtf="$download_dir/gencode.v${gencode}.primary_assembly.annotation.gtf"
processed_fa="$download_dir/sequence.fa"
processed_gtf="$download_dir/annotations.gtf"

function write_download_parameters() { 
    local gencode_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${gencode}/gencode.v${gencode}.primary_assembly.annotation.gtf.gz"
    local assembly_url="https://ftp.ensembl.org/pub/release-$ens_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

    echo "[\"$gencode_url\", \"$assembly_url\"]" > /tmp/output/genomegtfURLs
    mkdir -p $download_dir && echo "$download_dir" > /tmp/output/download_dir
}
function write_indices_parameters(){
  echo "$raw_fa" > /tmp/output/raw_fa
  echo "$raw_gtf" > /tmp/output/raw_gtf
  echo "$processed_fa" > /tmp/output/processed_fa
  echo "$processed_gtf" > /tmp/output/processed_gtf
  echo "$indices_dir" > /tmp/output/indices_dir
}

#check that one of aligners is chosen
[ -z "$useStar" ] && [ -z "$useKallisto" ] && [ -z "$useSalmon" ] && [ -z "$usePiscem" ] && echo "Must choose one aligner" && exit 1

if [ -n "$usePiscem" ] || [ -n "$useSalmon" ]; then
   [ -z "$useSplicei" ] && [ -z "$useSpliceu" ] && echo "Must choose one reference for salmon/piscem" && exit 1
fi

dir_keys=( "work_dir" "indices_dir" "download_dir" "aligneddir" "genome_dir" "counts_dir" )

for key in "${dir_keys[@]}"; do
  if [ -n "${!key}" ]; then
    echo "mkdir -p ${!key}"
  fi
done

write_download_parameters
write_indices_parameters

#declare -p
