#!/usr/bin/env python3
"""Test script for create_features_table.py"""

import os
import tempfile
from pathlib import Path

from create_features_table import create_features_table, read_barcode_list, read_assignments_file


def test_basic_functionality():
    """Test the basic functionality of create_features_table."""
    print("Testing basic functionality...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create test classification files
        test_files = {
            "ambiguous.txt": ["CELL_001", "CELL_002"],
            "doublets.txt": ["CELL_003", "CELL_001"],  # CELL_001 should be overwritten
            "missing_cells.txt": ["CELL_004"],
            "unassignable.txt": ["CELL_005", "CELL_006"]
        }
        
        for filename, barcodes in test_files.items():
            with open(temp_path / filename, 'w') as f:
                for barcode in barcodes:
                    f.write(f"{barcode}\n")
        
        # Create assignments.tsv file
        assignments_file = temp_path / "assignments.tsv"
        with open(assignments_file, 'w') as f:
            f.write("barcode\tvalue\n")  # Header
            f.write("CELL_007\t123\n")
            f.write("CELL_008\t456\n")
            f.write("CELL_001\t789\n")  # Should NOT overwrite existing classification
        
        # Output files
        output_file = temp_path / "features.tsv"
        categories_file = temp_path / "categories.txt"
        
        # Run the function
        create_features_table(
            str(temp_path),
            str(output_file),
            str(categories_file)
        )
        
        # Verify features.tsv output
        assert output_file.exists(), "features.tsv was not created"
        
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        # Check header
        assert lines[0].strip() == "barcode\tvalue", "Header is incorrect"
        
        # Parse content
        content = {}
        for line in lines[1:]:
            barcode, value = line.strip().split('\t')
            content[barcode] = value
        
        # Verify expected assignments
        expected = {
            "CELL_001": "doublet",      # Should be doublet (last classification file wins)
            "CELL_002": "ambiguous",
            "CELL_003": "doublet",
            "CELL_004": "missing_cells",
            "CELL_005": "unassignable",
            "CELL_006": "unassignable",
            "CELL_007": "BC123",        # From assignments.tsv
            "CELL_008": "BC456",        # From assignments.tsv
            # CELL_001 should NOT be BC789 because it was already classified
        }
        
        for barcode, expected_value in expected.items():
            assert barcode in content, f"Barcode {barcode} not found in output"
            assert content[barcode] == expected_value, \
                f"Expected {barcode}={expected_value}, got {content[barcode]}"
        
        # Verify categories file
        assert categories_file.exists(), "categories.txt was not created"
        
        with open(categories_file, 'r') as f:
            categories = [line.strip() for line in f.readlines()]
        
        # Check that it starts with fixed categories
        expected_start = ["ambiguous", "doublet", "missing_cells", "unassignable"]
        assert categories[:4] == expected_start, "Categories don't start with expected fixed categories"
        
        # Check that it includes BC categories
        assert "BC1" in categories, "BC1 not found in categories"
        assert "BC245979" in categories, "BC245979 not found in categories"
        assert len(categories) == 4 + 245979, f"Expected 245983 categories, got {len(categories)}"
        
        print("✓ Basic functionality test passed!")


def test_missing_files():
    """Test handling of missing files."""
    print("Testing missing files handling...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Only create one file
        with open(temp_path / "ambiguous.txt", 'w') as f:
            f.write("CELL_001\n")
        
        # No assignments.tsv file
        
        output_file = temp_path / "features.tsv"
        
        # Should not crash
        create_features_table(str(temp_path), str(output_file))
        
        # Verify output
        assert output_file.exists(), "features.tsv was not created"
        
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        assert len(lines) == 2, "Expected header + 1 data line"  # Header + CELL_001
        assert "CELL_001\tambiguous" in lines[1], "CELL_001 not classified correctly"
        
        print("✓ Missing files test passed!")


def test_file_reading_functions():
    """Test individual file reading functions."""
    print("Testing file reading functions...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Test barcode list reading
        barcode_file = temp_path / "test_barcodes.txt"
        with open(barcode_file, 'w') as f:
            f.write("CELL_001\n")
            f.write("CELL_002\n")
            f.write("# This is a comment\n")
            f.write("\n")  # Empty line
            f.write("CELL_003\n")
        
        barcodes = read_barcode_list(str(barcode_file))
        expected_barcodes = {"CELL_001", "CELL_002", "CELL_003"}
        assert barcodes == expected_barcodes, f"Expected {expected_barcodes}, got {barcodes}"
        
        # Test assignments reading
        assignments_file = temp_path / "test_assignments.tsv"
        with open(assignments_file, 'w') as f:
            f.write("barcode\tvalue\n")
            f.write("CELL_001\t123\n")
            f.write("CELL_002\t456\n")
            f.write("# Comment line\n")
            f.write("CELL_003\t789\n")
        
        assignments = read_assignments_file(str(assignments_file))
        expected_assignments = {
            "CELL_001": "BC123",
            "CELL_002": "BC456", 
            "CELL_003": "BC789"
        }
        assert assignments == expected_assignments, f"Expected {expected_assignments}, got {assignments}"
        
        print("✓ File reading functions test passed!")


def test_empty_files():
    """Test handling of empty files."""
    print("Testing empty files handling...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create empty files
        for filename in ["ambiguous.txt", "doublets.txt", "missing_cells.txt", "unassignable.txt"]:
            (temp_path / filename).touch()
        
        # Create empty assignments file with just header
        with open(temp_path / "assignments.tsv", 'w') as f:
            f.write("barcode\tvalue\n")
        
        output_file = temp_path / "features.tsv"
        
        # Should not crash
        create_features_table(str(temp_path), str(output_file))
        
        # Verify output (should just have header)
        assert output_file.exists(), "features.tsv was not created"
        
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        assert len(lines) == 1, "Expected only header line"
        assert lines[0].strip() == "barcode\tvalue", "Header is incorrect"
        
        print("✓ Empty files test passed!")


def run_all_tests():
    """Run all test functions."""
    print("Running tests for create_features_table.py")
    print("=" * 50)
    
    try:
        test_basic_functionality()
        test_missing_files()
        test_file_reading_functions()
        test_empty_files()
        
        print("=" * 50)
        print("🎉 All tests passed!")
        
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    run_all_tests()
