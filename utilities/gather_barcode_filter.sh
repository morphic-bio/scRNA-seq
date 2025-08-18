#!/bin/bash

# A script to organize alignment files based on gene directories.
# It takes five command-line arguments to specify input and output locations.

# --- Argument Validation ---
# Check if the correct number of arguments is provided.
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <gene_dir> <features_dir> <align_dir> <filter_file> <aligner>"
    echo "Example: $0 ./genes ./features ./alignments filtered_hits.txt bwa"
    exit 1
fi

# --- Assign Arguments to Variables ---
# Assign the command-line arguments to more descriptive variable names.
gene_dir="$1"
features_dir="$2"
align_dir="$3"
filter_file="$4"
aligner="$5"

# --- Directory Existence Check ---
# Ensure that the main source and destination directories exist before proceeding.
if [ ! -d "$gene_dir" ]; then
    echo "Error: Gene directory '$gene_dir' not found."
    exit 1
fi

if [ ! -d "$features_dir" ]; then
    echo "Error: Features directory '$features_dir' not found. Please create it first."
    exit 1
fi

if [ ! -d "$align_dir" ]; then
    echo "Error: Alignment directory '$align_dir' not found."
    exit 1
fi

echo "Starting file organization..."
echo "--------------------------------"

# --- Get Immediate Subdirectories from gene_dir ---
# Find all immediate subdirectories within gene_dir, get their base names,
# and store them in an array called 'subdirs'.
# The '*/' pattern ensures we only match directories.
subdirs=()
for d in "$gene_dir"/*/; do
    # Check if the path is a directory before adding it
    if [ -d "$d" ]; then
        subdirs+=("$(basename "$d")")
    fi
done

# --- Main Processing Loop ---
# Loop through each subdirectory name found.
for subdir_name in "${subdirs[@]}"; do
    echo "Processing subdirectory: $subdir_name"

    # Define the source file path
    source_file="$align_dir/$subdir_name/$aligner/$filter_file"

    # Define the destination directory path
    dest_dir="$features_dir/$subdir_name"

    # Create the destination subdirectory inside features_dir.
    # The -p flag ensures that the command doesn't fail if the directory already exists.
    mkdir -p "$dest_dir"
    echo "  -> Created or verified directory: $dest_dir"

    # Check if the source file exists before attempting to copy.
    if [ -f "$source_file" ]; then
        # Copy the specified filter file to the new subdirectory.
        cp "$source_file" "$dest_dir/"
        echo "  -> Copied '$source_file' to '$dest_dir/'"
    else
        echo "  -> Warning: Source file not found: $source_file"
    fi

    echo "--------------------------------"
done

echo "Script finished."
