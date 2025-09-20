#!/bin/bash

# Script to reassign larry barcodes to h5ad files using assign_larry_barcodes.py
# Usage: ./reassign_larry.sh [--dry-run]

set -euo pipefail

# Configuration
MSK_INPUT_DIR="/mnt/pikachu/MSK_30_KO"
ALIGNMENTS_DIR="/mnt/pikachu/storage/MSK-output-2/Alignments"
OUTPUT_DIR="/mnt/pikachu/new_MSK"
SCRIPT_DIR="/mnt/pikachu/scRNA-seq/utilities"
ASSIGN_SCRIPT="${SCRIPT_DIR}/assign_larry_barcodes.py"
NEW_COLUMN_NAME="larry_barcode"

# Parse command line arguments
DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "=== DRY RUN MODE - Commands will be echoed but not executed ==="
    echo
fi

# Check if the assign_larry_barcodes.py script exists
if [[ ! -f "$ASSIGN_SCRIPT" ]]; then
    echo "Error: assign_larry_barcodes.py not found at $ASSIGN_SCRIPT"
    exit 1
fi

# Check if input directory exists
if [[ ! -d "$MSK_INPUT_DIR" ]]; then
    echo "Error: Input directory not found: $MSK_INPUT_DIR"
    exit 1
fi

# Get first level subdirectories from MSK_30_KO and put just the names in an array
echo "Reading subdirectories from $MSK_INPUT_DIR..."
mapfile -t DIR_NAMES < <(find "$MSK_INPUT_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

# Check if any directories were found
if [[ ${#DIR_NAMES[@]} -eq 0 ]]; then
    echo "No subdirectories found in $MSK_INPUT_DIR"
    exit 1
fi

echo "Found ${#DIR_NAMES[@]} subdirectories:"
printf '  %s\n' "${DIR_NAMES[@]}"
echo

# Loop through each directory and process it
for dir_name in "${DIR_NAMES[@]}"; do
    echo "Processing directory: $dir_name"
    
    # Define input and output paths
    input_h5ad="${MSK_INPUT_DIR}/${dir_name}/counts.h5ad"
    mapping_file="${ALIGNMENTS_DIR}/${dir_name}/star/larry_features_em/features.tsv"
    categories_file="${ALIGNMENTS_DIR}/${dir_name}/star/larry_features_em/categories.txt"
    output_dir="${OUTPUT_DIR}/${dir_name}"
    output_h5ad="${output_dir}/counts.h5ad"
    
    # Check if required input files exist
    missing_files=()
    if [[ ! -f "$input_h5ad" ]]; then
        missing_files+=("$input_h5ad")
    fi
    if [[ ! -f "$mapping_file" ]]; then
        missing_files+=("$mapping_file")
    fi
    if [[ ! -f "$categories_file" ]]; then
        missing_files+=("$categories_file")
    fi
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        echo "  ⚠️  Skipping $dir_name - missing files:"
        printf '    %s\n' "${missing_files[@]}"
        echo
        continue
    fi
    
    # Build the command
    mkdir_cmd="mkdir -p \"$output_dir\""
    assign_cmd="python3 \"$ASSIGN_SCRIPT\" \"$input_h5ad\" \"$mapping_file\" \"$categories_file\" \"$NEW_COLUMN_NAME\" \"$output_h5ad\""
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY RUN] Would execute:"
        echo "    $mkdir_cmd"
        echo "    $assign_cmd"
    else
        echo "  Creating output directory: $output_dir"
        if mkdir -p "$output_dir"; then
            echo "  Executing: $assign_cmd"
            if python3 "$ASSIGN_SCRIPT" "$input_h5ad" "$mapping_file" "$categories_file" "$NEW_COLUMN_NAME" "$output_h5ad"; then
                echo "  ✓ Successfully processed $dir_name -> $output_h5ad"
            else
                echo "  ✗ Failed to process $dir_name" >&2
                # Continue with other directories instead of exiting
            fi
        else
            echo "  ✗ Failed to create output directory: $output_dir" >&2
        fi
    fi
    echo
done

if [[ "$DRY_RUN" == "true" ]]; then
    echo "=== DRY RUN COMPLETED - No files were actually created ==="
else
    echo "=== PROCESSING COMPLETED ==="
fi

echo "Summary: Processed ${#DIR_NAMES[@]} directories"
