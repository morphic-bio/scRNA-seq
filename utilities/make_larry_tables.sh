#!/bin/bash

# Script to find all larry_features_em directories and create features.tsv and categories.txt files
# Usage: ./make_larry_tables.sh [--dry-run]

set -euo pipefail

# Configuration
BASE_DIR="/mnt/pikachu/storage/MSK-output-2/Alignments"
SCRIPT_DIR="/mnt/pikachu/scRNA-seq/utilities"
CREATE_FEATURES_SCRIPT="${SCRIPT_DIR}/create_features_table.py"

# Parse command line arguments
DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "=== DRY RUN MODE - Commands will be echoed but not executed ==="
    echo
fi

# Check if the create_features_table.py script exists
if [[ ! -f "$CREATE_FEATURES_SCRIPT" ]]; then
    echo "Error: create_features_table.py not found at $CREATE_FEATURES_SCRIPT"
    exit 1
fi

# Find all larry_features_em directories and put them in an array
echo "Searching for larry_features_em directories in $BASE_DIR..."
mapfile -t LARRY_DIRS < <(find "$BASE_DIR" -type d -name "larry_features_em" | sort)

# Check if any directories were found
if [[ ${#LARRY_DIRS[@]} -eq 0 ]]; then
    echo "No larry_features_em directories found in $BASE_DIR"
    exit 1
fi

echo "Found ${#LARRY_DIRS[@]} larry_features_em directories:"
printf '  %s\n' "${LARRY_DIRS[@]}"
echo

# Loop through each directory and process it
for dir in "${LARRY_DIRS[@]}"; do
    echo "Processing directory: $dir"
    
    # Define output paths
    features_output="${dir}/features.tsv"
    categories_output="${dir}/categories.txt"
    
    # Build the command
    cmd="python3 \"$CREATE_FEATURES_SCRIPT\" \"$dir\" -o \"$features_output\" -c \"$categories_output\""
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY RUN] Would execute: $cmd"
    else
        echo "  Executing: $cmd"
        # Execute the command
        if python3 "$CREATE_FEATURES_SCRIPT" "$dir" -o "$features_output" -c "$categories_output"; then
            echo "  ✓ Successfully created features.tsv and categories.txt in $dir"
        else
            echo "  ✗ Failed to process $dir" >&2
            # Continue with other directories instead of exiting
        fi
    fi
    echo
done

if [[ "$DRY_RUN" == "true" ]]; then
    echo "=== DRY RUN COMPLETED - No files were actually created ==="
else
    echo "=== PROCESSING COMPLETED ==="
fi

echo "Summary: Processed ${#LARRY_DIRS[@]} directories"
