import anndata as ad
import os
import pandas as pd
import numpy as np

def _generate_ascii_plot(series, plot_type, max_width=40, max_categories=10, bins=20):
    """
    Generate an ASCII plot for a pandas Series.
    Internal helper function.
    """
    plot_str = []
    # Drop NaNs for plotting
    series = series.dropna()
    if series.empty:
        return ""

    if plot_type == 'histogram':
        # Use cut for uniform interval bins
        try:
            if series.nunique() < 2:
                return ""  # Not enough unique values for a histogram

            s_min, s_max = series.min(), series.max()

            # Define bin edges manually for better control
            # Use 0 as the start if data is non-negative, otherwise use the min
            start = 0 if s_min >= 0 else s_min
            
            # Use np.linspace to create uniform intervals ending exactly at max
            bin_edges = np.linspace(start, s_max, bins + 1)

            # Use pd.cut with manually defined bins
            # include_lowest=True ensures the first interval is left-inclusive
            labels = pd.cut(series, bins=bin_edges, include_lowest=True, duplicates='drop')
                
            counts = labels.value_counts().sort_index()

            if counts.empty:
                return ""

            max_val = counts.max()
            scale = max_width / max_val if max_val > 0 else 0

            plot_str.append("  " + "Histogram:")
            for interval, count in counts.items():
                bar = '█' * int(count * scale)
                # Format label to show correct interval closure, e.g., [0.00, 1.23] or (1.23, 2.46]
                left_bracket = '[' if interval.closed_left else '('
                right_bracket = ']' if interval.closed_right else ')'
                label = f"{left_bracket}{interval.left:.2f}, {interval.right:.2f}{right_bracket}"
                plot_str.append(f"  {label:>22} | {bar} ({count})")

        except Exception as e:
            return f"  Could not generate histogram: {str(e)[:100]}"

    elif plot_type == 'categorical_bar':
        counts = series.value_counts()
        
        title = "  Distribution:"
        if len(counts) > max_categories:
            title += f" (top {max_categories})"
            other_count = counts.iloc[max_categories:].sum()
            counts = counts.head(max_categories)
            if other_count > 0:
                # Use a dictionary to create a new series with '(other)'
                new_counts = counts.to_dict()
                new_counts['(other)'] = other_count
                counts = pd.Series(new_counts)
        
        plot_str.append(title)

        try:
            max_val = counts.max()
            scale = max_width / max_val if max_val > 0 else 0
            
            max_label_len = max(len(str(idx)) for idx in counts.index) if not counts.empty else 0

            for idx, count in counts.items():
                bar = '█' * int(count * scale)
                label = str(idx)
                plot_str.append(f"  {label:>{max_label_len}} | {bar} ({count})")
        except Exception as e:
            return f"  Could not generate bar plot: {str(e)[:100]}"

    return "\n".join(plot_str)

def infer_semantic_type(data, is_matrix=False, sample_size=1000):
    """
    Infer semantic type of data beyond just dtype.

    Parameters:
    -----------
    data : array-like or pandas Series
        Data to analyze
    is_matrix : bool
        Whether this is matrix data (for sampling)
    sample_size : int
        Number of values to sample for inference

    Returns:
    --------
    str
        Semantic type description
    """
    if is_matrix:
        # For matrices, sample some values
        if hasattr(data, 'nnz') and hasattr(data, 'data'):
            # Sparse matrix
            if data.nnz == 0:
                return "all zeros"
            
            # Use random sampling on the non-zero data to avoid bias from sequential slicing
            if data.nnz > sample_size:
                sample_indices = np.random.choice(data.nnz, sample_size, replace=False)
                sample_data = data.data[sample_indices]
            else:
                sample_data = data.data
        else:
            # Dense matrix (e.g., np.ndarray)
            flat_data = data.flatten()
            if len(flat_data) > sample_size:
                sample_indices = np.random.choice(len(flat_data), sample_size, replace=False)
                sample_data = flat_data[sample_indices]
            else:
                sample_data = flat_data
    else:
        # For series/1D data
        if len(data) > sample_size:
            sample_data = data.sample(sample_size) if hasattr(data, 'sample') else data[:sample_size]
        else:
            sample_data = data

    # Handle categorical data separately
    if hasattr(sample_data, 'dtype') and str(sample_data.dtype) == 'category':
        n_categories = len(sample_data.cat.categories)
        return f"categorical ({n_categories} categories)"

    # Handle object dtype (strings/mixed)
    if hasattr(sample_data, 'dtype') and sample_data.dtype == 'object':
        unique_values = sample_data.unique() if hasattr(sample_data, 'unique') else np.unique(sample_data)
        n_unique = len(unique_values)
        total_length = len(sample_data)

        if n_unique < total_length * 0.5:  # Less than 50% unique values
            return f"categorical ({n_unique} categories)"
        else:
            return "text/string data"

    # Remove NaN values for analysis
    if hasattr(sample_data, 'dropna'):
        clean_data = sample_data.dropna()
    else:
        # Convert to numpy array for consistent handling
        if hasattr(sample_data, 'values'):
            sample_array = sample_data.values
        else:
            sample_array = np.asarray(sample_data)

        if sample_array.dtype.kind in 'fc':  # float or complex
            clean_data = sample_array[~np.isnan(sample_array)]
        else:
            clean_data = sample_array

    if len(clean_data) == 0:
        return "all NaN/missing"

    # Convert to numpy array if it's still a pandas Series
    if hasattr(clean_data, 'values'):
        clean_data = clean_data.values

    unique_values, counts = np.unique(clean_data, return_counts=True)
    n_unique = len(unique_values)

    # Check for binary
    if n_unique == 2:
        # Floats with 2 unique values that are not {0, 1} are likely continuous, not binary.
        # This check defers to the more detailed float analysis that comes later.
        is_continuous_looking_float = (
            np.issubdtype(clean_data.dtype, np.floating) and
            set(unique_values) != {0.0, 1.0}
        )

        if not is_continuous_looking_float:
            val1, val2 = unique_values
            count1, count2 = counts

            # This must be checked first, as in Python {True, False} == {1, 0}
            if set(unique_values) == {True, False}:
                if is_matrix or isinstance(data, pd.Series):
                    try:
                        if hasattr(data, 'nnz'):  # Sparse matrix
                            count_true = np.sum(data.data)  # For boolean arrays, True is 1
                            count_false = (data.shape[0] * data.shape[1]) - data.nnz
                            return f"binary (True/False) (total counts: F≈{count_false}, T={count_true})"
                        else:  # Dense matrix or pd.Series
                            count_true = np.sum(data)
                            count_false = data.size - count_true
                            return f"binary (True/False) (total counts: F={count_false}, T={count_true})"
                    except Exception:
                        pass  # Fall through to sample counts on error
                # Sample-based fallback
                # np.unique sorts to [False, True], so counts[0] is for False.
                return f"binary (True/False) (counts in sample: F={counts[0]}, T={counts[1]})"
            elif set(unique_values) == {0, 1} or set(unique_values) == {0.0, 1.0}:
                if is_matrix or isinstance(data, pd.Series):
                    try:
                        if hasattr(data, 'nnz'):  # Sparse matrix
                            count_ones = np.sum(data.data == 1)
                            count_zeros = (data.shape[0] * data.shape[1]) - data.nnz
                            # Note: This assumes other non-zero values are negligible, as suggested by the sample.
                            return f"binary (0/1) (total counts: 0≈{count_zeros}, 1={count_ones})"
                        else:  # Dense matrix or pd.Series
                            count_ones = np.sum(data == 1)
                            count_zeros = np.sum(data == 0)
                            return f"binary (0/1) (total counts: 0={count_zeros}, 1={count_ones})"
                    except Exception:
                        pass  # Fall through to sample counts on error
                # Sample-based fallback for non-matrices or if full count fails
                return f"binary (0/1) (counts in sample: 0={count1}, 1={count2})"
            else:
                return f"binary ({val1}/{val2}) (counts in sample: {count1}/{count2})"

    # Now we can safely use numpy dtype checking
    try:
        # Check data type characteristics
        if np.issubdtype(clean_data.dtype, np.integer):
            if np.all(clean_data >= 0):
                if np.max(clean_data) <= 1:
                    # Catches samples of only 0s or only 1s.
                    if is_matrix or isinstance(data, pd.Series):
                        try:
                             if hasattr(data, 'nnz'):  # Sparse
                                count1 = np.sum(data.data == 1)
                                count0 = (data.shape[0] * data.shape[1]) - data.nnz
                                return f"binary/indicator integers (total counts: 0≈{count0}, 1={count1})"
                             else: # Dense or pd.Series
                                count1 = np.sum(data == 1)
                                count0 = np.sum(data == 0)
                                return f"binary/indicator integers (total counts: 0={count0}, 1={count1})"
                        except Exception:
                            pass # Fallback to sample counts
                    count1_sample = np.sum(clean_data)
                    count0_sample = len(clean_data) - count1_sample
                    return f"binary/indicator integers (counts in sample: 0={count0_sample}, 1={count1_sample})"
                else:
                    return "count data (non-negative integers)"
            else:
                return "signed integers"

        elif np.issubdtype(clean_data.dtype, np.floating):
            if np.all(clean_data >= 0):
                if np.all(clean_data <= 1):
                    return "normalized/probability values (0-1)"
                elif np.all(clean_data == np.round(clean_data)):
                    return "float-stored integers"
                else:
                    return "positive continuous values"
            else:
                # Check if might be log-transformed
                if np.all(clean_data > -20) and np.all(clean_data < 20):
                    return "continuous values (possibly log-transformed)"
                else:
                    return "continuous values"

        elif np.issubdtype(clean_data.dtype, np.bool_):
            if is_matrix or isinstance(data, pd.Series):
                try:
                    if hasattr(data, 'nnz'):  # Sparse
                        count_true = np.sum(data.data)
                        count_false = (data.shape[0] * data.shape[1]) - data.nnz
                        return f"boolean (total counts: F≈{count_false}, T={count_true})"
                    else: # Dense or pd.Series
                        count_true = np.sum(data)
                        count_false = data.size - count_true
                        return f"boolean (total counts: F={count_false}, T={count_true})"
                except Exception:
                    pass # Fallback to sample counts
            # Catches samples of only True or only False.
            count_true_sample = np.sum(clean_data)
            count_false_sample = len(clean_data) - count_true_sample
            return f"boolean (counts in sample: F={count_false_sample}, T={count_true_sample})"

        else:
            return f"other ({clean_data.dtype})"

    except Exception as e:
        # Fallback if dtype checking fails
        return f"unknown type (error: {str(e)[:50]})"

def analyze_data_types(data, name, is_matrix=False):
    """
    Analyze data types for a given dataset.

    Parameters:
    -----------
    data : array-like or pandas DataFrame
        Data to analyze
    name : str
        Name of the dataset
    is_matrix : bool
        Whether this is matrix data

    Returns:
    --------
    str
        Type analysis string
    """
    if is_matrix:
        base_dtype = str(data.dtype)
        semantic_type = infer_semantic_type(data, is_matrix=True)
        return f"{base_dtype} | {semantic_type}"
    else:
        # For DataFrame columns
        return infer_semantic_type(data, is_matrix=False)

def summarize_h5ad(file_path, show_samples=True, n_samples=5, analyze_types=True, show_plots=True):
    """
    Read an AnnData h5ad file and summarize its structure.

    Parameters:
    -----------
    file_path : str
        Path to the h5ad file
    show_samples : bool
        Whether to show sample values from obs and main matrix
    n_samples : int
        Number of sample values to display (default: 5)
    analyze_types : bool
        Whether to perform detailed type analysis
    show_plots : bool
        Whether to generate ASCII plots for obs/var annotations

    Returns:
    --------
    str
        A text report summarizing the structure of the AnnData object
    """
    # Check if file exists
    if not os.path.exists(file_path):
        return f"Error: File '{file_path}' not found."

    # Load the AnnData object
    try:
        adata = ad.read_h5ad(file_path)
    except Exception as e:
        return f"Error loading file: {str(e)}"

    # Build report
    report = []
    report.append(f"AnnData Summary: {os.path.basename(file_path)}")
    report.append("=" * 50)

    # Basic information
    report.append(f"Dimensions: {adata.shape[0]} cells × {adata.shape[1]} genes/features")

    # Main matrix X info
    x_type = type(adata.X).__name__
    x_dtype = str(adata.X.dtype)
    report.append(f"\nMain Matrix (X): {x_type} of shape {adata.shape}")
    report.append(f"Data type: {x_dtype}")

    if analyze_types:
        try:
            x_semantic = infer_semantic_type(adata.X, is_matrix=True)
            report.append(f"Inferred type: {x_semantic}")
        except Exception as e:
            report.append(f"Type inference error: {str(e)[:100]}")

    # Check for sparsity
    if hasattr(adata.X, 'data') and hasattr(adata.X, 'indices'):
        sparsity = 1.0 - (len(adata.X.data) / (adata.shape[0] * adata.shape[1]))
        report.append(f"Sparsity: {sparsity:.2%}")

    # Show sample values from main matrix
    if show_samples:
        report.append(f"\nSample values from main matrix (first {min(n_samples, adata.shape[0])} cells, first {min(n_samples, adata.shape[1])} genes):")

        # Get sample subset
        n_cells_sample = min(n_samples, adata.shape[0])
        n_genes_sample = min(n_samples, adata.shape[1])

        # Convert to dense if sparse for sampling
        try:
            if hasattr(adata.X, 'toarray'):
                sample_data = adata.X[:n_cells_sample, :n_genes_sample].toarray()
            else:
                sample_data = adata.X[:n_cells_sample, :n_genes_sample]

            # Create a mini dataframe for display
            sample_df = pd.DataFrame(sample_data, 
                                   index=[f"Cell_{i}" for i in range(n_cells_sample)],
                                   columns=[f"Gene_{i}" for i in range(n_genes_sample)])

            # Format the output nicely
            report.append(sample_df.to_string(float_format='%.3f'))
        except Exception as e:
            report.append(f"Error displaying sample values: {str(e)}")

    # Layers
    if adata.layers:
        report.append("\nLayers:")
        for layer_name, layer in adata.layers.items():
            layer_info = f"- {layer_name}: {type(layer).__name__} of shape {layer.shape}, dtype: {layer.dtype}"
            
            # Add min/max values
            try:
                # Handle all-zero sparse matrices, where .min()/.max() can fail
                if hasattr(layer, 'nnz') and layer.nnz == 0:
                    min_val, max_val = 0.0, 0.0
                else:
                    min_val = layer.min()
                    max_val = layer.max()
                layer_info += f", range: [{min_val:.3f}, {max_val:.3f}]"
            except Exception:
                # Fallback for types that don't support min/max, or empty arrays
                pass

            if analyze_types:
                try:
                    layer_semantic = infer_semantic_type(layer, is_matrix=True)
                    layer_info += f", inferred: {layer_semantic}"
                except Exception as e:
                    layer_info += f", type inference error: {str(e)[:50]}"
            report.append(layer_info)

            if show_samples:
                try:
                    report.append(f"  Sample values from {layer_name}:")
                    n_cells_sample = min(n_samples, layer.shape[0])
                    n_genes_sample = min(3, layer.shape[1])  # Show fewer genes for layers

                    if hasattr(layer, 'toarray'):
                        sample_layer = layer[:n_cells_sample, :n_genes_sample].toarray()
                    else:
                        sample_layer = layer[:n_cells_sample, :n_genes_sample]

                    sample_layer_df = pd.DataFrame(sample_layer,
                                                 index=[f"Cell_{i}" for i in range(n_cells_sample)],
                                                 columns=[f"Gene_{i}" for i in range(n_genes_sample)])
                    report.append("  " + sample_layer_df.to_string(float_format='%.3f').replace('\n', '\n  '))
                except Exception as e:
                    report.append(f"  Error displaying layer samples: {str(e)}")

    # Observations (obs)
    if adata.obs.columns.size > 0:
        report.append("\nObservation Annotations (obs):")
        for col in adata.obs.columns:
            try:
                n_unique = adata.obs[col].nunique()
                base_type = str(adata.obs[col].dtype)
                semantic_type_full = ""
                type_info = f"- {col}: {base_type}, {n_unique} unique values"

                if analyze_types:
                    semantic_type_full = infer_semantic_type(adata.obs[col], is_matrix=False)
                    type_info += f" | {semantic_type_full}"
                report.append(type_info)

                if show_plots and semantic_type_full:
                    semantic_type_short = semantic_type_full.split('(')[0].strip()
                    plot_str = ""
                    if 'binary' in semantic_type_short or 'boolean' in semantic_type_short or semantic_type_short == 'categorical':
                        plot_str = _generate_ascii_plot(adata.obs[col], 'categorical_bar')
                    elif adata.obs[col].nunique() > 1 and pd.api.types.is_numeric_dtype(adata.obs[col].dtype):
                        plot_str = _generate_ascii_plot(adata.obs[col], 'histogram')
                    
                    if plot_str:
                        report.append(plot_str)

                if show_samples:
                    sample_values = adata.obs[col].head(n_samples)
                    if adata.obs[col].dtype == 'object' or str(adata.obs[col].dtype) == 'category':
                        values_str = ", ".join([f"'{v}'" for v in sample_values.astype(str)])
                    else:
                        values_str = ", ".join([f"{v}" for v in sample_values])
                    report.append(f"  Sample values: {values_str}")
            except Exception as e:
                report.append(f"- {col}: Error analyzing column: {str(e)}")

    # Variables (var)
    if adata.var.columns.size > 0:
        report.append("\nVariable Annotations (var):")
        for col in adata.var.columns:
            try:
                n_unique = adata.var[col].nunique()
                base_type = str(adata.var[col].dtype)
                semantic_type_full = ""
                type_info = f"- {col}: {base_type}, {n_unique} unique values"

                if analyze_types:
                    semantic_type_full = infer_semantic_type(adata.var[col], is_matrix=False)
                    type_info += f" | {semantic_type_full}"
                report.append(type_info)

                if show_plots and semantic_type_full:
                    semantic_type_short = semantic_type_full.split('(')[0].strip()
                    plot_str = ""
                    if 'binary' in semantic_type_short or 'boolean' in semantic_type_short or semantic_type_short == 'categorical':
                        plot_str = _generate_ascii_plot(adata.var[col], 'categorical_bar')
                    elif adata.var[col].nunique() > 1 and pd.api.types.is_numeric_dtype(adata.var[col].dtype):
                        plot_str = _generate_ascii_plot(adata.var[col], 'histogram')

                    if plot_str:
                        report.append(plot_str)

                if show_samples:
                    sample_values = adata.var[col].head(n_samples)
                    if adata.var[col].dtype == 'object' or str(adata.var[col].dtype) == 'category':
                        values_str = ", ".join([f"'{v}'" for v in sample_values.astype(str)])
                    else:
                        values_str = ", ".join([f"{v}" for v in sample_values])
                    report.append(f"  Sample values: {values_str}")
            except Exception as e:
                report.append(f"- {col}: Error analyzing column: {str(e)}")

    # Remaining components
    for name, component in [
        ("Unstructured Annotations (uns)", adata.uns),
        ("Observation Matrices (obsm)", adata.obsm),
        ("Variable Matrices (varm)", adata.varm),
        ("Observation Pairwise Matrices (obsp)", adata.obsp),
        ("Variable Pairwise Matrices (varp)", adata.varp)
    ]:
        if component:
            report.append(f"\n{name}:")
            for key, value in component.items():
                try:
                    shape_info = f"shape {value.shape}" if hasattr(value, "shape") else ""
                    value_type = type(value).__name__

                    type_info = f"- {key}: {value_type} {shape_info}"
                    if hasattr(value, "dtype"):
                        type_info += f", dtype: {value.dtype}"
                        if analyze_types and hasattr(value, "shape"):
                            semantic_type = infer_semantic_type(value, is_matrix=True)
                            type_info += f" | {semantic_type}"

                    report.append(type_info)

                    # Show samples for matrices
                    if show_samples and hasattr(value, "shape") and len(value.shape) >= 1:
                        if len(value.shape) == 1:  # 1D array
                            sample_vals = value[:min(n_samples, len(value))]
                            vals_str = ", ".join([f"{v:.3f}" if isinstance(v, (int, float, np.number)) else str(v) for v in sample_vals])
                            report.append(f"  Sample values: {vals_str}")
                        elif len(value.shape) == 2:  # 2D array
                            n_rows = min(3, value.shape[0])
                            n_cols = min(3, value.shape[1])
                            sample_2d = value[:n_rows, :n_cols]
                            report.append(f"  Sample values ({n_rows}×{n_cols} subset):")
                            sample_2d_df = pd.DataFrame(sample_2d)
                            report.append("  " + sample_2d_df.to_string(float_format='%.3f').replace('\n', '\n  '))
                except Exception as e:
                    report.append(f"- {key}: Error analyzing component: {str(e)}")

    # ASCII representation of structure
    report.append("\nASCII Structure:")
    report.append("AnnData")
    report.append("├── X: main data matrix")

    if adata.layers:
        report.append("├── layers")
        last_layer = list(adata.layers.keys())[-1]
        for layer in adata.layers.keys():
            prefix = "│   └──" if layer == last_layer else "│   ├──"
            report.append(f"{prefix} {layer}")

    report.append("├── obs: cell annotations")
    if adata.obs.columns.size > 0:
        for i, col in enumerate(adata.obs.columns):
            prefix = "│   └──" if i == len(adata.obs.columns) - 1 else "│   ├──"
            report.append(f"{prefix} {col}")

    report.append("├── var: gene annotations")
    if adata.var.columns.size > 0:
        for i, col in enumerate(adata.var.columns):
            prefix = "│   └──" if i == len(adata.var.columns) - 1 else "│   ├──"
            report.append(f"{prefix} {col}")

    # Add remaining elements to ASCII tree
    remaining = []
    if adata.uns: remaining.append("uns")
    if adata.obsm: remaining.append("obsm")
    if adata.varm: remaining.append("varm")
    if adata.obsp: remaining.append("obsp")
    if adata.varp: remaining.append("varp")

    for i, elem in enumerate(remaining):
        is_last = i == len(remaining) - 1
        prefix = "└──" if is_last else "├──"
        parent_prefix = "    " if is_last else "│   "

        component_map = {
            "uns": (adata.uns, "uns: unstructured annotations"),
            "obsm": (adata.obsm, "obsm: observation matrices"),
            "varm": (adata.varm, "varm: variable matrices"),
            "obsp": (adata.obsp, "obsp: observation pairwise matrices"),
            "varp": (adata.varp, "varp: variable pairwise matrices")
        }

        component, label = component_map[elem]
        if component:
            report.append(f"{prefix} {label}")
            keys = list(component.keys())
            for j, key in enumerate(keys):
                key_prefix = f"{parent_prefix}└──" if j == len(keys) - 1 else f"{parent_prefix}├──"
                report.append(f"{key_prefix} {key}")

    return "\n".join(report)

# Example usage with command line arguments
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Summarize AnnData h5ad file structure')
    parser.add_argument('file_path', help='Path to the h5ad file')
    parser.add_argument('--no-samples', action='store_true', 
                       help='Disable showing sample values')
    parser.add_argument('--no-types', action='store_true',
                       help='Disable detailed type analysis')
    parser.add_argument('--no-plots', action='store_true',
                       help='Disable generating ASCII plots for annotations')
    parser.add_argument('--n-samples', type=int, default=5,
                       help='Number of sample values to show (default: 5)')

    args = parser.parse_args()

    summary = summarize_h5ad(args.file_path, 
                           show_samples=not args.no_samples, 
                           n_samples=args.n_samples,
                           analyze_types=not args.no_types,
                           show_plots=not args.no_plots)
    print(summary)
