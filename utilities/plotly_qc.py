import anndata as ad
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import argparse

def plot_interactive_qc_histograms(adata: ad.AnnData, output_html: str = "qc_histograms.html"):
    """
    Creates an interactive HTML file with histograms for numeric columns in adata.obs.

    A dropdown menu allows the user to select which column to display.
    Buttons allow toggling between all values and non-zero values.

    Parameters:
    -----------
    adata : anndata.AnnData
        The AnnData object containing the data.
    output_html : str
        The path to save the output HTML file.
    """
    # Find all numeric columns in the observation metadata
    numeric_cols = adata.obs.select_dtypes(include=np.number).columns.tolist()

    if not numeric_cols:
        print("No numeric columns found in adata.obs to plot.")
        return

    # Create a Plotly figure
    fig = go.Figure()

    # Pre-calculate both all-value and non-zero data lists to make updates fast
    all_data = [adata.obs[col] for col in numeric_cols]
    nonzero_data = [d[d > 0] for d in all_data]

    # Add a histogram trace for each numeric column, starting with the full data
    for i, col in enumerate(numeric_cols):
        fig.add_trace(go.Histogram(
            x=all_data[i],
            name=col,
            # Make the first trace visible by default
            visible=(i == 0)
        ))

    # --- Create Menus ---

    # 1. Dropdown for variable selection
    dropdown_buttons = []
    for i, col in enumerate(numeric_cols):
        # Create a visibility list: True for the current trace, False for all others
        visibility = [False] * len(numeric_cols)
        visibility[i] = True
        
        dropdown_buttons.append(dict(
            label=col,
            method="update",
            args=[
                {"visible": visibility},
                {
                    "title_text": f"Distribution of {col}",
                    "xaxis.title.text": col,
                    "yaxis.autorange": True
                }
            ]
        ))

    # 2. Buttons to toggle between all and non-zero values
    scale_buttons = [
        dict(
            label="All Values",
            method="restyle",
            args=[{"x": all_data}]
        ),
        dict(
            label="Non-Zero Values",
            method="restyle",
            args=[{"x": nonzero_data}]
        )
    ]

    # Update the figure layout with both menus
    fig.update_layout(
        updatemenus=[
            # Dropdown menu for variable
            dict(
                active=0,
                buttons=dropdown_buttons,
                direction="down",
                pad={"r": 10, "t": 10},
                x=0.05, xanchor="left", y=1.15, yanchor="top"
            ),
            # Buttons for data scale
            dict(
                type="buttons",
                direction="right",
                active=0,
                buttons=scale_buttons,
                pad={"r": 10, "t": 10},
                x=0.45, xanchor="left", y=1.15, yanchor="top"
            )
        ],
        title_text=f"Distribution of {numeric_cols[0]}",
        xaxis_title=numeric_cols[0],
        yaxis_title="Count",
        bargap=0.1,
        annotations=[
            dict(text="QC Metric:", showarrow=False,
                 x=0, y=1.13, yref="paper", xref="paper", align="left"),
            dict(text="Data:", showarrow=False,
                 x=0.4, y=1.13, yref="paper", xref="paper", align="left")
        ]
    )

    # Save the figure to an HTML file
    fig.write_html(output_html)
    print(f"Interactive plot saved to {output_html}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate interactive QC plots from an AnnData file.')
    parser.add_argument('file_path', help='Path to the h5ad file')
    parser.add_argument('--output', '-o', default='qc_histograms.html',
                       help='Output HTML file name (default: qc_histograms.html)')
    
    args = parser.parse_args()

    try:
        adata = ad.read_h5ad(args.file_path)
        plot_interactive_qc_histograms(adata, args.output)
    except Exception as e:
        print(f"Error: Could not process the file. {e}")
