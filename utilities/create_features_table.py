#!/usr/bin/env python3
"""
Create features table from barcode classification files.

This script processes barcode classification files and creates a features.tsv file
with barcode-value pairs, along with an optional categories file.
"""

import argparse
import os
from pathlib import Path
from typing import Dict, Set


def read_barcode_list(file_path: str) -> Set[str]:
    """Read a file containing one barcode per line and return a set of barcodes.
    
    Args:
        file_path: Path to the file containing barcodes
        
    Returns:
        Set of barcode strings (stripped of whitespace)
    """
    barcodes = set()
    if os.path.exists(file_path):
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                barcode = line.strip()
                if barcode and not barcode.startswith('#'):  # Skip empty lines and comments
                    barcodes.add(barcode)
        print(f"Read {len(barcodes)} barcodes from {file_path}")
    else:
        print(f"File not found: {file_path}")
    return barcodes


def read_assignments_file(file_path: str) -> Dict[str, str]:
    """Read assignments.tsv file with barcode and value columns.
    
    Args:
        file_path: Path to the assignments.tsv file
        
    Returns:
        Dictionary mapping barcode to BC-prefixed value
    """
    assignments = {}
    if os.path.exists(file_path):
        with open(file_path, 'r', encoding='utf-8') as f:
            first_line = True
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Skip header line if it looks like a header
                if first_line:
                    first_line = False
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        # Check if this looks like a header
                        barcode_col = parts[0].strip().lower()
                        value_col = parts[1].strip().lower()
                        header_tokens = {'barcode', 'barcodes', 'cell', 'cell_id', 'cellid', 
                                       'value', 'values', 'assignment', 'label', 'category'}
                        if barcode_col in header_tokens or value_col in header_tokens:
                            continue  # Skip header line
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    barcode = parts[0].strip()
                    value = parts[1].strip()
                    if barcode:
                        assignments[barcode] = f"BC{value}"
                else:
                    print(f"Warning: Line {line_num} in {file_path} does not have 2 tab-separated columns: {line}")
        
        print(f"Read {len(assignments)} assignments from {file_path}")
    else:
        print(f"File not found: {file_path}")
    
    return assignments


def create_features_table(input_directory: str, output_file: str = "features.tsv", 
                         categories_file: str = None) -> None:
    """Create features table from barcode classification files.
    
    Args:
        input_directory: Directory containing the classification files
        output_file: Output file path for features.tsv
        categories_file: Optional categories file path
    """
    input_dir = Path(input_directory)
    
    # Initialize barcode dictionary
    barcode_dict = {}
    
    # Define classification files in order (later files can overwrite earlier ones)
    classification_files = [
        ("ambiguous.txt", "ambiguous"),
        ("doublets.txt", "doublet"),
        ("missing_cells.txt", "missing_cells"),
        ("unassignable.txt", "unassignable")
    ]
    
    # Read classification files
    print("Reading classification files...")
    for filename, label in classification_files:
        file_path = input_dir / filename
        barcodes = read_barcode_list(str(file_path))
        
        # Add barcodes to dictionary (overwriting if they already exist)
        for barcode in barcodes:
            barcode_dict[barcode] = label
    
    print(f"Total barcodes after classification files: {len(barcode_dict)}")
    
    # Read assignments file
    assignments_file = input_dir / "assignments.tsv"
    assignments = read_assignments_file(str(assignments_file))
    
    # Add assignments to dictionary (but don't overwrite existing values)
    initial_count = len(barcode_dict)
    for barcode, bc_value in assignments.items():
        if barcode not in barcode_dict:
            barcode_dict[barcode] = bc_value
    
    print(f"Added {len(barcode_dict) - initial_count} new barcodes from assignments.tsv")
    print(f"Total barcodes in final dictionary: {len(barcode_dict)}")
    
    # Write features.tsv file
    print(f"Writing features table to {output_file}...")
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("barcode\tvalue\n")  # Header
        for barcode in sorted(barcode_dict.keys()):
            f.write(f"{barcode}\t{barcode_dict[barcode]}\n")
    
    print(f"Wrote {len(barcode_dict)} barcode-value pairs to {output_file}")
    
    # Create categories file if requested
    if categories_file:
        create_categories_file(categories_file, barcode_dict)


def create_categories_file(categories_file: str, barcode_dict: Dict[str, str]) -> None:
    """Create categories file with all possible category values.
    
    Args:
        categories_file: Path to the categories file to create
        barcode_dict: Dictionary of barcode assignments to extract BC values from
    """
    print(f"Creating categories file: {categories_file}")
    
    # Start with fixed categories
    categories = ["ambiguous", "doublet", "missing_cells", "unassignable"]
    
    # Add BC categories from 1 to 245979
    bc_categories = [f"BC{i}" for i in range(1, 245980)]  # 1 to 245979 inclusive
    categories.extend(bc_categories)
    
    # Write categories file
    with open(categories_file, 'w', encoding='utf-8') as f:
        for category in categories:
            f.write(f"{category}\n")
    
    print(f"Wrote {len(categories)} categories to {categories_file}")
    
    # Report which categories are actually used
    used_categories = set(barcode_dict.values())
    unused_count = len(categories) - len(used_categories)
    print(f"Categories actually used in data: {len(used_categories)}")
    print(f"Unused categories: {unused_count}")


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Create features table from barcode classification files"
    )
    parser.add_argument(
        "input_directory",
        help="Directory containing classification files (ambiguous.txt, doublets.txt, etc.)"
    )
    parser.add_argument(
        "-o", "--output",
        default="features.tsv",
        help="Output file path for features table (default: features.tsv)"
    )
    parser.add_argument(
        "-c", "--categories",
        help="Create categories file at specified path"
    )
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.isdir(args.input_directory):
        print(f"Error: Input directory does not exist: {args.input_directory}")
        return 1
    
    try:
        create_features_table(
            input_directory=args.input_directory,
            output_file=args.output,
            categories_file=args.categories
        )
        print("✓ Features table creation completed successfully!")
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
