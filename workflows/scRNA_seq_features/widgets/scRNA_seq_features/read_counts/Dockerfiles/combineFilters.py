#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import anndata as ad
import os
from scipy.sparse import issparse
from scipy.sparse import csr_matrix
import scanpy as sc
import scipy.sparse as sp

def add_filter_to_anndata(adata, barcode_file, filter_column_name, present_value=True):
    """
    Updates the AnnData object in place with a binary filter column in obs based on the presence of barcodes.

    Parameters:
    adata (AnnData): The input AnnData object.
    barcode_file (str): Path to the text file containing barcodes (one per line).
    filter_column_name (str): The name of the binary filter column to be added to obs.
    present_value (bool): The value (True/False) to assign if the barcode is present.
    """
    print(f"\n--- Starting add_filter_to_anndata ---")
    print(f"Barcode file: {barcode_file}")
    print(f"Filter column name: {filter_column_name}")
    print(f"Present value: {present_value}")

    # Read barcodes into a set
    with open(barcode_file, 'r') as f:
        barcodes = {line.strip() for line in f}

    # Fast, hash-based membership test
    matches = adata.obs_names.isin(barcodes)  # returns a boolean array

    # If present_value=True we want 1 for matches, else 0 for matches
    if present_value:
        adata.obs[filter_column_name] = matches.astype(int)
    else:
        adata.obs[filter_column_name] = (~matches).astype(int)

    print(f"--- Finished add_filter_to_anndata for {filter_column_name} ---")

def calculate_adaptive_thresholds(adata, singlet_mask, mt_genes, n_mad=3):
    """
    Calculate adaptive thresholds based on median and MAD (Median Absolute Deviation).
    Now only calculates the max_genes threshold adaptively.

    Parameters:
    adata (AnnData): The input AnnData object.
    singlet_mask (pd.Series): Boolean mask for singlet cells.
    mt_genes (list): List of mitochondrial genes.
    n_mad (float): Number of MADs to use for threshold calculation.

    Returns:
    int: max_genes threshold
    """
    print(f"\n--- Starting calculate_adaptive_thresholds ---")
    print(f"Using n_mad: {n_mad}")
    print(f"Number of cells in input adata: {adata.n_obs}")
    print(f"Number of singlet cells for calculation: {np.sum(singlet_mask)}")

    # Calculate metrics for singlet cells only
    print("Subsetting AnnData to singlet cells for adaptive threshold calculation...")
    singlet_mask_bool = singlet_mask.astype(bool)
    singlet_adata = adata[singlet_mask_bool]

    # Calculate number of genes per cell
    print("Calculating number of genes per cell for singlets...")
    n_genes = np.asarray(np.sum(singlet_adata.X > 0, axis=1)).flatten()
    print(f"Calculated n_genes for {len(n_genes)} cells.")

    # Calculate median and MAD for genes
    print("Calculating median and MAD for n_genes...")
    n_genes_median = np.median(n_genes)
    n_genes_mad = np.median(np.abs(n_genes - n_genes_median))
    print(f"n_genes_median: {n_genes_median}, n_genes_mad: {n_genes_mad}")

    # Calculate max_genes threshold
    print("Calculating adaptive threshold for max_genes only...")
    max_genes = int(n_genes_median + n_mad * n_genes_mad)
    
    print(f"\nAdaptive threshold calculation:")
    print(f"Number of genes - Median: {n_genes_median:.1f}, MAD: {n_genes_mad:.1f}")
    print(f"Calculated max_genes threshold: {max_genes}")

    print(f"--- Finished calculate_adaptive_thresholds ---")
    return max_genes

def filter_cells(adata, min_genes, max_genes, mt_genes, mt_pct_cutoff, filter_name='filter'):
    """
    Filters cells based on the number of genes and percentage of mitochondrial genes.

    Parameters:
    adata (AnnData): The input AnnData object.
    min_genes (int): Minimum number of genes a cell must have to be kept.
    max_genes (int): Maximum number of genes a cell can have to be kept.
    mt_genes (list): List of mitochondrial genes.
    mt_pct_cutoff (float): Maximum allowed percentage of mitochondrial genes.
    filter_name (str): The name of the binary filter column to be added to obs.

    Returns:
    AnnData: The AnnData object with the binary filter column added to obs.
    """
    print(f"\n--- Starting filter_cells ---")
    print(f"Applying filters with parameters:")
    print(f"  min_genes: {min_genes}")
    print(f"  max_genes: {max_genes}")
    print(f"  mt_pct_cutoff: {mt_pct_cutoff}")
    print(f"  filter_name: {filter_name}")
    print(f"  Number of mitochondrial genes provided: {len(mt_genes)}")

    # Calculate the total number of genes and mitochondrial genes for each cell
    print("Calculating 'n_genes' per cell...")
    adata.obs['n_genes'] = np.asarray(np.sum(adata.X > 0, axis=1)).flatten()
    # Calculate the total counts of mitochondrial genes for each cell
    print("Calculating 'mt_counts' and 'total_counts' per cell...")
    if issparse(adata.X):
        print("Input data is sparse. Performing memory-efficient calculations for mt_pct in filter_cells.")
        adata.obs['mt_counts'] = np.asarray(adata[:, mt_genes].X.sum(axis=1)).flatten()
        adata.obs['total_counts'] = np.asarray(adata.X.sum(axis=1)).flatten()
    else:
        print("Input data is dense for mt_pct calculation in filter_cells.")
        adata.obs['mt_counts'] = np.sum(adata[:, mt_genes].X, axis=1)
        adata.obs['total_counts'] = np.sum(adata.X, axis=1)

    # Calculate the percentage of mitochondrial genes for each cell
    print("Calculating 'mt_pct' per cell...")
    # Avoid division by zero: replace 0s in total_counts with 1s for calculation,
    # then explicitly set mt_pct to 0 where total_counts was 0.
    total_counts_series = adata.obs['total_counts']
    mt_counts_series = adata.obs['mt_counts']
    
    # Create a safe denominator (replace 0 with 1 to avoid division by zero error)
    total_counts_safe = total_counts_series.copy()
    total_counts_safe[total_counts_safe == 0] = 1
    
    adata.obs['mt_pct'] = (mt_counts_series / total_counts_safe) * 100
    # Where original total_counts was 0, mt_pct should be 0
    adata.obs.loc[total_counts_series == 0, 'mt_pct'] = 0

    # Apply the filter conditions
    print("Applying filter conditions...")
    filter_array = (adata.obs['n_genes'] >= min_genes) & (adata.obs['n_genes'] <= max_genes) & (adata.obs['mt_pct'] <= mt_pct_cutoff)

    print(f"Adding '{filter_name}' column to AnnData.obs...")
    adata.obs[filter_name] = filter_array
    print(f"Number of cells passing this filter: {np.sum(filter_array)}")
    print(f"--- Finished filter_cells ---")

def add_doublet_scores_to_anndata(adata, doublet_score_file):
    print(f"\n--- Starting add_doublet_scores_to_anndata ---")
    print(f"Reading doublet score file: {doublet_score_file}")
    doublet_scores_df = pd.read_csv(doublet_score_file, sep='\t')
    print(f"Read {len(doublet_scores_df)} scores from file.")
    print("Setting 'Barcode' as index for doublet scores DataFrame...")
    doublet_scores_df.set_index('Barcode', inplace=True)
    print("Mapping doublet scores to AnnData.obs['doublet_scores']...")
    adata.obs['doublet_scores'] = adata.obs_names.map(doublet_scores_df['Score'])
    print(f"Number of scores mapped: {adata.obs['doublet_scores'].notna().sum()}")
    print(f"--- Finished add_doublet_scores_to_anndata ---")

def create_qc_plot(adata, directory, min_genes, max_genes, mt_pct_cutoff, num_singlet_filtered, num_single_cells, n_mad=None):
    """
    Create QC plot with gene distribution histogram and MT percentage.
    
    Parameters:
    adata (AnnData): The AnnData object
    directory (str): Output directory for saving plots
    min_genes (int): Minimum genes threshold
    max_genes (int): Maximum genes threshold
    mt_pct_cutoff (float): MT percentage cutoff
    num_singlet_filtered (int): Number of cells passing all filters
    num_single_cells (int): Total number of singlet cells
    n_mad (float): Number of MADs used for adaptive filtering (optional)
    """
    print("Creating interactive gene distribution histogram with quantiles...")
    import plotly.graph_objects as go
    import numpy as np
    from plotly.subplots import make_subplots
    
    singlet_cells = adata.obs['singlet']
    singlet_data = adata[singlet_cells]
    
    # Also get the filtered data
    filtered_cells = adata.obs['singlet_filtered']
    filtered_data = adata[filtered_cells]
    
    # Calculate quantiles for binning (20 quantiles) based on singlet data
    gene_counts = singlet_data.obs['n_genes']
    mt_pcts = singlet_data.obs['mt_pct']
    filtered_gene_counts = filtered_data.obs['n_genes']
    
    # Start quantiles from 0 instead of minimum
    max_gene_count = gene_counts.max()
    quantiles = np.linspace(0, max_gene_count, 21)
    
    # Create the figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    # Pre-calculate bin heights to find the max height for y-axis scaling
    bin_counts, bin_edges = np.histogram(
        gene_counts, 
        bins=20, 
        range=(0, max_gene_count)
    )
    max_bin_height = bin_counts.max()
    
    # Add histogram trace for all singlet cells
    hist1 = go.Histogram(
        x=gene_counts,
        xbins=dict(
            start=0,
            end=max_gene_count,
            size=max_gene_count/20
        ),
        marker=dict(
            color='blue',
            line=dict(color='white', width=1)
        ),
        opacity=0.7,
        name='All singlets'
    )
    fig.add_trace(hist1, secondary_y=False)
    
    # Add histogram trace for filtered cells
    hist2 = go.Histogram(
        x=filtered_gene_counts,
        xbins=dict(
            start=0,
            end=max_gene_count,
            size=max_gene_count/20
        ),
        marker=dict(
            color='green',
            line=dict(color='white', width=1)
        ),
        opacity=0.7,
        name='Filtered singlets'
    )
    fig.add_trace(hist2, secondary_y=False)
    
    # Calculate mean MT percentage per gene count quantile bin
    mt_per_bin = []
    bin_centers = []
    
    # For each bin, calculate the mean mt_pct
    for i in range(len(bin_edges)-1):
        # Get the cells that fall within this bin
        mask = (gene_counts >= bin_edges[i]) & (gene_counts < bin_edges[i+1])
        if np.sum(mask) > 0:
            mt_per_bin.append(mt_pcts[mask].mean())
        else:
            mt_per_bin.append(0)
        bin_centers.append((bin_edges[i] + bin_edges[i+1]) / 2)
    
    # Add MT percentage line on secondary axis
    fig.add_trace(
        go.Scatter(
            x=bin_centers,
            y=mt_per_bin,
            mode='lines+markers',
            marker=dict(
                color='red',
                size=8
            ),
            line=dict(
                color='red',
                width=2
            ),
            name='MT % (mean)'
        ),
        secondary_y=True
    )
    
    # Add vertical lines for thresholds - different colors for each
    fig.add_trace(
        go.Scatter(
            x=[min_genes, min_genes],
            y=[0, max_bin_height * 1.1],
            mode='lines',
            line=dict(color='red', width=2, dash='dash'),
            name=f'Min: {min_genes}'
        ),
        secondary_y=False
    )
    
    fig.add_trace(
        go.Scatter(
            x=[max_genes, max_genes],
            y=[0, max_bin_height * 1.1],
            mode='lines',
            line=dict(color='orange', width=2, dash='dash'),
            name=f'Max: {max_genes}'
        ),
        secondary_y=False
    )
    
    # Add horizontal line for MT percentage cutoff
    fig.add_trace(
        go.Scatter(
            x=[0, max_gene_count],
            y=[mt_pct_cutoff, mt_pct_cutoff],
            mode='lines',
            line=dict(color='purple', width=2, dash='dash'),
            name=f'MT% cutoff: {mt_pct_cutoff}%'
        ),
        secondary_y=True
    )

    # Create title with MAD information if available
    title_text = "Gene Distribution and MT Percentage by Quantile"
    if n_mad is not None:
        title_text += f" (Max genes: {n_mad} MADs)"

    # Update layout for better visualization - centered title
    fig.update_layout(
        height=600,
        width=900,
        title=dict(
            text=title_text,
            x=0.5,  # Center the title
            xanchor='center'
        ),
        xaxis_title="Number of genes",
        margin=dict(l=50, r=50, t=80, b=50),  # Increased top margin from 50 to 80
        showlegend=True,
        template='plotly_white',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        barmode='overlay'  # Overlay histograms instead of stacking
    )
    
    # Update y-axes titles
    fig.update_yaxes(title_text="Number of cells", range=[0, max_bin_height * 1.1], secondary_y=False)
    fig.update_yaxes(title_text="Mitochondrial percentage (%)", range=[0, max(mt_pcts.max() * 1.1, mt_pct_cutoff * 1.2)], secondary_y=True)
    
    # Add annotation with thresholds and filtering stats - include MAD info if available
    stats_text = f"Min: {min_genes} | Max: {max_genes} | Passed: {num_singlet_filtered}/{num_single_cells}"
    if n_mad is not None:
        stats_text += f" | MADs: {n_mad}"
    
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.01, y=0.99,
        text=stats_text,
        showarrow=False,
        font=dict(size=10),
        align="left",
        bgcolor="white",
        bordercolor="black",
        borderwidth=1
    )

    # Calculate and display min/max/median values
    median_genes = np.median(gene_counts)
    filtered_median_genes = np.median(filtered_gene_counts)
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.01, y=0.93,
        text=f"All: median={int(median_genes)} | Filtered: median={int(filtered_median_genes)}",
        showarrow=False,
        font=dict(size=10),
        align="left",
        bgcolor="white",
        bordercolor="black",
        borderwidth=1
    )

    print("Saving plot files...")
    output_html = os.path.join(directory, 'gene_quantile_histogram.html')
    output_png = os.path.join(directory, 'gene_quantile_histogram.png')
    
    fig.write_html(output_html)
    fig.write_image(output_png, scale=2)  # Higher resolution for the static image
    
    print(f"Plots saved to: {output_html} and {output_png}")

def ensure_csr_format(adata):
    """Ensure AnnData X matrix and all layers are in CSR format."""
    # Fix main matrix
    if hasattr(adata.X, 'toarray') and not sp.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()
    
    # Fix layers if they exist
    if hasattr(adata, 'layers'):
        for layer_name in adata.layers:
            layer_matrix = adata.layers[layer_name]
            if hasattr(layer_matrix, 'toarray') and not sp.isspmatrix_csr(layer_matrix):
                adata.layers[layer_name] = layer_matrix.tocsr()
    
    return adata

#parse arguments to get input file
print("Parsing command line arguments...")
parser = argparse.ArgumentParser(description='Create a filtered AnnData object')
parser.add_argument('--input_file', type=str, help='Input AnnData file', required=True)
#add flag for mintochondrial genes file
parser.add_argument('--mito_genes', type=str, help='File containing mitochondrial genes')
#add flag for adaptive filtering
parser.add_argument('--adaptive_filter', action='store_true', help='Use adaptive filtering based on median and MAD')
#add flag for number of MADs
parser.add_argument('--n_mad', type=float, default=3, help='Number of MADs to use for adaptive filtering (default: 3)')

args = parser.parse_args()
print(f"Arguments parsed: {args}")

#if the mito genes file is provided read it in (one per line)
if args.mito_genes:
    print(f"Reading mitochondrial genes from file: {args.mito_genes}")
    with open(args.mito_genes, 'r') as f:
        mito_genes = {line.strip() for line in f}
    print(f"Read {len(mito_genes)} mitochondrial genes.")
else:
    mito_genes = ['ENSG00000198888','ENSG00000198763','ENSG00000198804','ENSG00000198712','ENSG00000228253','ENSG00000198899','ENSG00000198938','ENSG00000198840','ENSG00000212907','ENSG00000198886','ENSG00000198786','ENSG00000198695','ENSG00000198727']
    print(f"Using default list of {len(mito_genes)} mitochondrial genes.")

input_file = args.input_file
directory = os.path.dirname(input_file)
print(f"Input AnnData file: {input_file}")
print(f"Output directory: {directory}")

print(f"\nReading AnnData object from {input_file}...")
adata=ad.read_h5ad(input_file)
print(f"Successfully read AnnData object with {adata.n_obs} observations and {adata.n_vars} variables.")

print("\nAdding filter for non-empty barcodes...")
add_filter_to_anndata(adata, os.path.join(directory, 'non_empty_barcodes.txt'), 'non_empty', present_value=True)
print("\nAdding filter for doublet barcodes...")
add_filter_to_anndata(adata, os.path.join(directory, 'doublet_barcodes.txt'), 'doublet', present_value=True)
print("\nAdding doublet scores...")
add_doublet_scores_to_anndata(adata, os.path.join(directory, 'filtered_barcodes_with_scores.txt'))

# Calculate singlet cells
print("\nCalculating singlet cells (non_empty AND NOT doublet)...")
singlet_cells = adata.obs['non_empty'] & ~adata.obs['doublet'] 
adata.obs['singlet'] = singlet_cells.astype(bool)
print(f"Number of singlet cells identified: {np.sum(singlet_cells)}")

# Determine filtering thresholds
print("\nDetermining filtering thresholds...")
if args.adaptive_filter:
    print("\nUsing adaptive filtering for max_genes only...")
    # Get min_genes and mt_pct_cutoff from environment or use defaults
    min_genes = os.environ.get('min_genes', 200)
    min_genes = int(min_genes)
    
    mt_pct_cutoff = os.environ.get('mt_pct_cutoff', 5.0)
    mt_pct_cutoff = float(mt_pct_cutoff)
    
    # Calculate max_genes adaptively
    max_genes = calculate_adaptive_thresholds(
        adata, 
        singlet_cells, 
        list(mito_genes),
        n_mad=args.n_mad
    )
    
    print(f"\nUsing thresholds:")
    print(f"  Min genes: {min_genes} (from input/default)")
    print(f"  Max genes: {max_genes} (adaptive)")
    print(f"  MT % cutoff: {mt_pct_cutoff}% (from input/default)")
else:
    # Use fixed thresholds from environment variables
    print("Attempting to use fixed thresholds from environment variables...")
    min_genes = os.environ.get('min_genes', 200)
    max_genes = os.environ.get('max_genes', 2500)
    mt_pct_cutoff = os.environ.get('mt_pct_cutoff', 5.0)  # Get from env with default 5%
    print(f"Raw values from env (or default): min_genes='{min_genes}', max_genes='{max_genes}', mt_pct_cutoff='{mt_pct_cutoff}'")
    #convert values to appropriate types
    min_genes = int(min_genes)
    max_genes = int(max_genes)
    mt_pct_cutoff = float(mt_pct_cutoff)
    print(f"\nUsing fixed thresholds:")
    print(f"  Min genes: {min_genes}")
    print(f"  Max genes: {max_genes}")
    print(f"  MT % cutoff: {mt_pct_cutoff}%")

# Apply filtering
print("\nApplying cell quality filters...")
filter_cells(adata, min_genes, max_genes, list(mito_genes), mt_pct_cutoff, 'filter') # Ensure mito_genes is a list

# Calculate combined filter
print("\nCalculating combined filter (singlet AND passed QC filter)...")
singlet_filtered=adata.obs['singlet'] & adata.obs['filter']
adata.obs['singlet_filtered'] = singlet_filtered
print(f"Number of cells in combined 'singlet_filtered' group: {np.sum(singlet_filtered)}")

# Print statistics
num_non_empty = np.sum(adata.obs['non_empty'])
num_single_cells = np.sum(singlet_cells)
num_filtered = np.sum(adata.obs['filter'])
num_singlet_filtered = np.sum(singlet_filtered)

print(f"\nCell filtering statistics:")
print(f"Number of non-empty cells: {num_non_empty}")
print(f"Number of single cells: {num_single_cells}")
print(f"Number of cells that passed the filter: {num_filtered}")
print(f"Number of single cells that passed the filter: {num_singlet_filtered}")

# Save QC plots if using adaptive filtering
if args.adaptive_filter:
    print("\nAdaptive filtering was used, generating simplified quantile-based histogram...")
    create_qc_plot(adata, directory, min_genes, max_genes, mt_pct_cutoff, num_singlet_filtered, num_single_cells, args.n_mad)
else:
    print("\nAdaptive filtering not used, skipping QC plot generation.")

print (adata)   

output_file = os.path.join(directory, 'unfiltered_counts.h5ad')
print(f"\nWriting raw AnnData object (with all .obs columns) to {output_file}...")
# Ensure CSR format before writing to prevent AnnData warnings
adata = ensure_csr_format(adata)
adata.write(output_file)
print(f"Successfully wrote {output_file}")

#create a binary filter from adata.obs['singlet_filtered'] and use it to subset adata
print("\nCreating filtered AnnData object based on 'singlet_filtered' mask...")
adata_filtered = adata[adata.obs['singlet_filtered']].copy() # Corrected subsetting
print(f"Filtered AnnData object created with {adata_filtered.n_obs} observations.")
output_file = os.path.join(directory, 'filtered_counts.h5ad')
print(f"Writing filtered AnnData object to {output_file}...")
# Ensure CSR format before writing to prevent AnnData warnings
adata_filtered = ensure_csr_format(adata_filtered)
adata_filtered.write(output_file)
print(f"Successfully wrote {output_file}")
print(adata_filtered)
print("\n--- Script finished ---")
