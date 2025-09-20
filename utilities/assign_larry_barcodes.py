import argparse
import os
import re
from typing import List

import anndata as ad
import pandas as pd


def _read_categories_file(categories_file_path: str) -> List[str]:
    """Read a categories file into an ordered list of category strings.

    Supports either one category per line, or a single line with comma/tab/space-separated values.
    Lines beginning with '#' and empty lines are ignored.
    """
    with open(categories_file_path, "r", encoding="utf-8") as handle:
        raw = handle.read()

    # Remove comments and trim
    lines = []
    for line in raw.splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            lines.append(line)

    # If a single line contains separators, split it
    if len(lines) == 1 and re.search(r"[\t, ]", lines[0]):
        split_values = re.split(r"[\t, ]+", lines[0])
        categories = [v for v in (s.strip() for s in split_values) if v]
    else:
        categories = lines

    # Deduplicate while preserving order
    seen = set()
    ordered_unique: List[str] = []
    for c in categories:
        if c not in seen:
            seen.add(c)
            ordered_unique.append(c)

    if not ordered_unique:
        raise ValueError(f"No categories found in file: {categories_file_path}")

    return ordered_unique


def _read_two_column_mapping_file(mapping_file_path: str) -> pd.Series:
    """Read a two-column mapping file (barcode, label) into a Series indexed by barcode.

    - Accepts TSV/CSV/whitespace delimited files automatically (sep=None, engine='python').
    - Handles optional header row gracefully; if the first row looks like a header, it is skipped.
    - Keeps the last occurrence for duplicate barcodes.
    - Returns a pd.Series with index=barcode (str), values=label (str or NA).
    """
    # Read without assuming header to avoid dropping first data row when header is absent
    df = pd.read_csv(
        mapping_file_path,
        sep=None,
        engine="python",
        dtype=str,
        header=None,
        comment="#",
    )

    if df.shape[1] < 2:
        raise ValueError(
            f"Expected at least 2 columns in mapping file, found {df.shape[1]}: {mapping_file_path}"
        )

    # For whitespace-separated files, pandas may create extra columns filled with NaN
    # Take the first column (barcode) and the last non-empty column (label)
    barcode_col = df.iloc[:, 0]
    
    # Find the last column that has non-null values for most rows
    label_col = None
    for col_idx in range(df.shape[1] - 1, 0, -1):  # Start from last column, go backwards
        col_data = df.iloc[:, col_idx]
        non_null_count = col_data.notna().sum()
        if non_null_count > 0:  # Found a column with actual data
            label_col = col_data
            break
    
    if label_col is None:
        # Fallback to second column if no non-null column found
        label_col = df.iloc[:, 1]
    
    df = pd.DataFrame({"barcode": barcode_col, "label": label_col})

    # Detect and drop header row if present
    header_barcode_tokens = {
        "barcode",
        "barcodes",
        "cell",
        "cell_id",
        "cellid",
        "cellbarcode",
        "cell_barcode",
        "cb",
        "index",
        "obs",
        "obs_names",
    }
    header_label_tokens = {
        "label",
        "labels",
        "category",
        "categories",
        "group",
        "groups",
        "class",
        "type",
        "cluster",
        "assignment",
    }
    if len(df) > 0:
        first_barcode = str(df.iloc[0, 0]).strip().lower()
        first_label = str(df.iloc[0, 1]).strip().lower()
        if first_barcode in header_barcode_tokens and (
            first_label in header_label_tokens or first_label in header_barcode_tokens
        ):
            df = df.iloc[1:].copy()

    # Normalize whitespace
    df["barcode"] = df["barcode"].astype(str).str.strip()
    df["label"] = df["label"].astype(str).str.strip()

    # Drop rows with empty barcodes
    df = df[df["barcode"].astype(bool)]

    if df.empty:
        raise ValueError(f"No valid barcode-label rows found in mapping file: {mapping_file_path}")

    # Keep the last occurrence for duplicate barcodes
    df = df.drop_duplicates(subset=["barcode"], keep="last")

    series = df.set_index("barcode")["label"]
    return series


def assign_larry_barcodes(
    input_h5ad_path: str,
    mapping_file_path: str,
    categories_file_path: str,
    new_obs_column_name: str,
    output_h5ad_path: str,
) -> None:
    """Assign categorical labels from a two-column mapping file to an AnnData object's obs.

    This function:
      1) Loads the input .h5ad AnnData
      2) Creates a new obs column (named by new_obs_column_name) initialized to NA
      3) Intersects barcodes between the mapping file and AnnData.obs_names
      4) Fills matched barcodes with the corresponding label
      5) Ensures the obs column is a categorical dtype with categories from the categories file
      6) Writes the updated AnnData to output_h5ad_path
    """
    if not os.path.exists(input_h5ad_path):
        raise FileNotFoundError(f"Input h5ad not found: {input_h5ad_path}")
    if not os.path.exists(mapping_file_path):
        raise FileNotFoundError(f"Mapping file not found: {mapping_file_path}")
    if not os.path.exists(categories_file_path):
        raise FileNotFoundError(f"Categories file not found: {categories_file_path}")

    # Load AnnData
    adata = ad.read_h5ad(input_h5ad_path)

    # Read categories and mapping
    categories = _read_categories_file(categories_file_path)
    mapping_series = _read_two_column_mapping_file(mapping_file_path)

    # Prepare a new Series of NA values for all cells
    obs_index = pd.Index(adata.obs_names.astype(str))
    result_series = pd.Series(pd.NA, index=obs_index, dtype="object")

    # Intersect barcodes and assign
    mapping_index = pd.Index(mapping_series.index.astype(str))
    intersecting_barcodes = obs_index.intersection(mapping_index)
    if len(intersecting_barcodes) == 0:
        # Still create the column as all NA with the provided categories
        print(
            f"No barcodes in mapping file intersect AnnData.obs_names. Creating NA column '{new_obs_column_name}'."
        )
    else:
        result_series.loc[intersecting_barcodes] = mapping_series.loc[intersecting_barcodes].values

    # Cast to categorical with provided categories (values not in categories become NA)
    categorical_dtype = pd.api.types.CategoricalDtype(categories=categories, ordered=False)
    result_series = result_series.astype(categorical_dtype)

    # Assign to obs (overwrite if exists)
    if new_obs_column_name in adata.obs.columns:
        print(f"Column '{new_obs_column_name}' exists in obs and will be overwritten.")
    adata.obs[new_obs_column_name] = result_series

    # Write to new file
    adata.write_h5ad(output_h5ad_path)

    print(
        f"Wrote updated AnnData with column '{new_obs_column_name}' to: {output_h5ad_path}. "
        f"Matched {len(intersecting_barcodes)} of {adata.n_obs} cells."
    )


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Assign categorical labels from a two-column mapping file to a new obs column in an AnnData .h5ad."
        )
    )
    parser.add_argument("input_h5ad", help="Path to input .h5ad AnnData file")
    parser.add_argument(
        "mapping_file",
        help=(
            "Two-column file with barcodes in first column and categorical label in second column. "
            "Supports CSV/TSV/whitespace. Optional header is handled."
        ),
    )
    parser.add_argument(
        "categories_file",
        help=(
            "File containing list of categories (one per line, or a single line of comma/tab/space-separated categories)."
        ),
    )
    parser.add_argument("new_column_name", help="Name of the new obs column to create")
    parser.add_argument("output_h5ad", help="Path to write the updated .h5ad AnnData file")
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    assign_larry_barcodes(
        input_h5ad_path=args.input_h5ad,
        mapping_file_path=args.mapping_file,
        categories_file_path=args.categories_file,
        new_obs_column_name=args.new_column_name,
        output_h5ad_path=args.output_h5ad,
    )


if __name__ == "__main__":
    main()


