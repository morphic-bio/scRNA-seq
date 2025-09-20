#!/usr/bin/env python3
"""Test script for assign_larry_barcodes.py"""

import os
import tempfile
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from assign_larry_barcodes import assign_larry_barcodes, _read_categories_file, _read_two_column_mapping_file


def create_test_anndata(barcodes, n_genes=100):
    """Create a test AnnData object with specified barcodes."""
    n_cells = len(barcodes)
    
    # Create random expression matrix
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))
    
    # Create obs and var dataframes
    obs = pd.DataFrame(index=barcodes)
    obs['existing_column'] = np.random.choice(['A', 'B', 'C'], n_cells)
    
    var = pd.DataFrame(index=[f'Gene_{i}' for i in range(n_genes)])
    var['gene_type'] = 'protein_coding'
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata


def test_basic_functionality():
    """Test the basic functionality of assign_larry_barcodes."""
    print("Testing basic functionality...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create test data
        test_barcodes = ['CELL_001', 'CELL_002', 'CELL_003', 'CELL_004', 'CELL_005']
        adata = create_test_anndata(test_barcodes)
        
        # Save test h5ad
        input_h5ad = temp_path / "test_input.h5ad"
        adata.write_h5ad(input_h5ad)
        
        # Create mapping file (some barcodes match, some don't)
        mapping_file = temp_path / "mapping.tsv"
        mapping_data = [
            "barcode\tlabel",  # header
            "CELL_001\tType1",
            "CELL_003\tType2", 
            "CELL_999\tType1",  # This barcode doesn't exist in anndata
            "CELL_005\tType3"
        ]
        with open(mapping_file, 'w') as f:
            f.write('\n'.join(mapping_data))
        
        # Create categories file
        categories_file = temp_path / "categories.txt"
        with open(categories_file, 'w') as f:
            f.write("Type1\nType2\nType3\nType4")  # Type4 won't be used
        
        # Output file
        output_h5ad = temp_path / "test_output.h5ad"
        
        # Run the function
        assign_larry_barcodes(
            str(input_h5ad),
            str(mapping_file),
            str(categories_file),
            "cell_type",
            str(output_h5ad)
        )
        
        # Verify results
        result_adata = ad.read_h5ad(output_h5ad)
        
        # Check that new column exists
        assert "cell_type" in result_adata.obs.columns, "New column not created"
        
        # Check that it's categorical
        assert isinstance(result_adata.obs["cell_type"].dtype, pd.CategoricalDtype), "Column is not categorical"
        
        # Check categories
        expected_categories = ['Type1', 'Type2', 'Type3', 'Type4']
        actual_categories = list(result_adata.obs["cell_type"].cat.categories)
        assert actual_categories == expected_categories, f"Categories mismatch: {actual_categories} vs {expected_categories}"
        
        # Check specific assignments
        assert result_adata.obs.loc["CELL_001", "cell_type"] == "Type1"
        assert result_adata.obs.loc["CELL_003", "cell_type"] == "Type2"
        assert result_adata.obs.loc["CELL_005", "cell_type"] == "Type3"
        assert pd.isna(result_adata.obs.loc["CELL_002", "cell_type"])  # Should be NA
        assert pd.isna(result_adata.obs.loc["CELL_004", "cell_type"])  # Should be NA
        
        print("✓ Basic functionality test passed!")


def test_no_intersection():
    """Test case where no barcodes intersect."""
    print("Testing no intersection case...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create test data
        test_barcodes = ['CELL_001', 'CELL_002', 'CELL_003']
        adata = create_test_anndata(test_barcodes)
        
        # Save test h5ad
        input_h5ad = temp_path / "test_input.h5ad"
        adata.write_h5ad(input_h5ad)
        
        # Create mapping file with non-matching barcodes
        mapping_file = temp_path / "mapping.tsv"
        with open(mapping_file, 'w') as f:
            f.write("CELL_999\tType1\n")
            f.write("CELL_998\tType2\n")
        
        # Create categories file
        categories_file = temp_path / "categories.txt"
        with open(categories_file, 'w') as f:
            f.write("Type1\nType2")
        
        # Output file
        output_h5ad = temp_path / "test_output.h5ad"
        
        # Run the function
        assign_larry_barcodes(
            str(input_h5ad),
            str(mapping_file),
            str(categories_file),
            "cell_type",
            str(output_h5ad)
        )
        
        # Verify results - all should be NA
        result_adata = ad.read_h5ad(output_h5ad)
        assert all(pd.isna(result_adata.obs["cell_type"])), "Expected all values to be NA"
        
        print("✓ No intersection test passed!")


def test_categories_file_formats():
    """Test different formats for categories files."""
    print("Testing categories file formats...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Test comma-separated categories
        categories_file = temp_path / "categories_comma.txt"
        with open(categories_file, 'w') as f:
            f.write("Type1, Type2, Type3")
        
        categories = _read_categories_file(str(categories_file))
        assert categories == ['Type1', 'Type2', 'Type3'], f"Comma-separated parsing failed: {categories}"
        
        # Test tab-separated categories
        categories_file = temp_path / "categories_tab.txt"
        with open(categories_file, 'w') as f:
            f.write("Type1\tType2\tType3")
        
        categories = _read_categories_file(str(categories_file))
        assert categories == ['Type1', 'Type2', 'Type3'], f"Tab-separated parsing failed: {categories}"
        
        # Test one-per-line categories with comments
        categories_file = temp_path / "categories_lines.txt"
        with open(categories_file, 'w') as f:
            f.write("# This is a comment\nType1\nType2  # inline comment\n\nType3\n")
        
        categories = _read_categories_file(str(categories_file))
        assert categories == ['Type1', 'Type2', 'Type3'], f"Line-by-line parsing failed: {categories}"
        
        print("✓ Categories file format tests passed!")


def test_mapping_file_formats():
    """Test different formats for mapping files."""
    print("Testing mapping file formats...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Test CSV format with header
        mapping_file = temp_path / "mapping.csv"
        with open(mapping_file, 'w') as f:
            f.write("barcode,label\n")
            f.write("CELL_001,Type1\n")
            f.write("CELL_002,Type2\n")
        
        mapping = _read_two_column_mapping_file(str(mapping_file))
        assert mapping.loc['CELL_001'] == 'Type1', f"Expected Type1, got {mapping.loc['CELL_001']}"
        assert mapping.loc['CELL_002'] == 'Type2', f"Expected Type2, got {mapping.loc['CELL_002']}"
        
        # Test TSV format without header
        mapping_file = temp_path / "mapping.tsv"
        with open(mapping_file, 'w') as f:
            f.write("CELL_001\tType1\n")
            f.write("CELL_002\tType2\n")
        
        mapping = _read_two_column_mapping_file(str(mapping_file))
        assert mapping.loc['CELL_001'] == 'Type1', f"Expected Type1, got {mapping.loc['CELL_001']}"
        assert mapping.loc['CELL_002'] == 'Type2', f"Expected Type2, got {mapping.loc['CELL_002']}"
        
        # Test whitespace-separated format
        mapping_file = temp_path / "mapping.txt"
        with open(mapping_file, 'w') as f:
            f.write("CELL_001   Type1\n")
            f.write("CELL_002   Type2\n")
        
        mapping = _read_two_column_mapping_file(str(mapping_file))
        assert mapping.loc['CELL_001'] == 'Type1', f"Expected Type1, got {mapping.loc['CELL_001']}"
        assert mapping.loc['CELL_002'] == 'Type2', f"Expected Type2, got {mapping.loc['CELL_002']}"
        
        print("✓ Mapping file format tests passed!")


def test_duplicate_barcodes():
    """Test handling of duplicate barcodes in mapping file."""
    print("Testing duplicate barcode handling...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create mapping file with duplicates (should keep last)
        mapping_file = temp_path / "mapping.tsv"
        with open(mapping_file, 'w') as f:
            f.write("CELL_001\tType1\n")
            f.write("CELL_001\tType2\n")  # This should be kept
            f.write("CELL_002\tType3\n")
        
        mapping = _read_two_column_mapping_file(str(mapping_file))
        assert mapping.loc['CELL_001'] == 'Type2', "Should keep last occurrence of duplicate"
        assert mapping.loc['CELL_002'] == 'Type3'
        
        print("✓ Duplicate barcode test passed!")


def run_all_tests():
    """Run all test functions."""
    print("Running tests for assign_larry_barcodes.py")
    print("=" * 50)
    
    try:
        test_basic_functionality()
        test_no_intersection()
        test_categories_file_formats()
        test_mapping_file_formats()
        test_duplicate_barcodes()
        
        print("=" * 50)
        print("🎉 All tests passed!")
        
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        raise


if __name__ == "__main__":
    run_all_tests()
