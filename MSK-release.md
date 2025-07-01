# Release notes
## TL;DR
1. Added cell types (user and predicted)
2. Changed to MADS (median + standard deviation) to find outliers for number of genes to filter
3. Bugfix for MT%
4. The top occuring feature is now shown with the counts of the top two features
5. Star is the only aligner now
6. New QC visualizations
7. Gencode 44 to match with other experiments
8. Feature counts from gene expression files now added as a layer 

### Organization of files
The files are organized by experiments (A,B,C,D,E,F,G_1,G_2,H,I,J,L_1,L_2) with features.h5ad counts.h5ad and QC directory

### Structure of files
counts.h5ad
```
├── X: main data matrix
├── layers
│   ├── ambiguous
│   ├── denoised
│   ├── spliced
│   └── unspliced
├── obs: cell annotations
│   ├── is_cell
│   ├── non_empty
│   ├── doublet
│   ├── doublet_scores
│   ├── singlet
│   ├── n_genes
│   ├── mt_counts
│   ├── total_counts
│   ├── mt_pct
│   ├── filter
│   ├── singlet_filtered
│   ├── is_featured
│   ├── best_feature
│   ├── feature1_count
│   ├── feature2_count
│   ├── is_denoised
│   ├── user_celltype
│   └── predicted_celltype
├── var: gene annotations
```
features.h5ad
```
├── X: main data matrix
├── layers
│   ├── features
│   └── gex_features
├── obs: cell annotations
│   ├── is_expressed
│   ├── is_cell
│   ├── non_empty
│   ├── doublet
│   ├── doublet_scores
│   ├── singlet
│   ├── n_genes
│   ├── mt_counts
│   ├── total_counts
│   ├── mt_pct
│   ├── filter
│   ├── singlet_filtered
│   └── is_featured
├── var: gene annotations

```

### Definition of X, obs, var, layers
#### X (Data Matrix): (<number_of_barcodes>, 36601)
This is the basic 2D matrix with gene counts or feature counts. The gene counts include the feature counts as they are expressed. These are at the columns that follow the normal gene columns. The number of rows depends on the aligner. The gene counts include intron counts (as per CellRanger >= 7.0) and is the sum of the ambiguous, spliced and unspliced layers.
#### obs: (<number_of_barcodes>, 12)
##### is_cell
Cellranger and Star-solo run their own empty-cell detection. This has a value of 1 when the internal methodology identifies the barcode to be a non-empty cell.
##### non_empty
The EmptyDropsCellRanger [algorithm](https://github.com/MarioniLab/DropletUtils/pull/64/commits/5a0b6c2f91b4767b92e3fb0894a2c05c9eaa22d4) is run to detect empty cells. A value of 1 indicates that the cell is non-empty
##### doublet
[scDblFinder](https://github.com/plger/scDblFinder) was run to detect doublets. A value of 1 indicates that the barcode is a doublet
##### doublet_score
Confidence that the cell is a doublet. 
#### n_genes
Number of genes with non-zero counts
##### mt_counts
Total number of counts from mitochondrial genes
##### mt_pct
Percentage of counts that are from mitochondrial genes
##### filter
Has value of 1 if mt percentage is < 5 and the number of genes is < median + 3 standard deviations
##### singlet
Has a value of 1 if the barcode is not detected as a doublet and it is not empty
##### singlet_filtered
Has a value of 1 if the barcode is a singlet and passes the filter
#### var: (38592, 0) (for gene counts file)
##### unspliced, spliced, ambiguous
These are the different count values (layers) for the barcode, gene count matrix. The unspliced, spliced and ambiguous counts are defined by this figure from the kallisto preprint  https://www.biorxiv.org/content/10.1101/2022.12.02.518832v2.full with M being mature/spliced, A, ambiguous and N, nascent/unspliced. The definitions are those defined by Velocyto for all the aligners except Kallisto. Note that the workaround is to assign Ambiguous counts to spliced counts, in which case all aligners have the same definition.
##### best_feature
The feature with the most counts - empty if there is a tie
##### feature1_count
The number of counts for the best feature
##### feature2_count
The number of counts for the second best feature
##### denoised
The denoised layer is produced by cellBender a deep-learning package from Broad that detects and removes ambient counts. Note that cellBender uses an initial empty-cell detection method and only calculates the denoised counts for non-empty cells. The rest of the layer is zero-filled. Note that it is possible but unlikely that some non-empty cells detected by the EmptyDrops algorithm may produce zero counts when denoised. This is unlikely as cellBender intentionally uses a very generous initial criterion for filtering out empty cells.
#### For feature counts file
####  X (Data Matrix): (<number_of_barcodes, 105)
This is the basic 2D matrix with feature counts. The counts are determined using the Biodepot software from the targeted sequencing. UMI deduping proceeds by finding the feature that has the most reads assigned to a barcode-umi. Matches of up to 3 ambiguous bases and a max Hamming distance of 5 were considered as long as the match was unambiguous. Sequence barcodes are fuzzy matched using the CellRanger methodology. 
##### expression
The feature barcodes (larry barcodes) are expressed and the scRNA-seq fastqs can be searched for matches directly. The same matching methods were used except that de-duping was done using the standard RNA-seq strategy - all barcode-umi counts that mapped to a feature are reduced to 1.
##### sc_layer
This uses the alignment to the GFP-feature sequence to determine the counts. Previously, the feature sequence alone was used which gave reasonable counts only with Kallisto. The more complete sequence should work with star and piscem though this has not been tested.

### Feature counts
The main layer is from the targeted sequencing. The gex layer is from using the Process_features tool to find the features from the expression fastq files.

### QC folders
These are at the top level and inside the sample directories. The html files are interactive plotly graphs.

#### Feature QC
#### Run stats
In the output directory some statistics from the run are kept in "stats.txt"
```
Total feature counts 31993392
Total deduped feature counts 6765455
Total unique barcode UMIs 9148050
Total whitelisted barcodes 259694
Total_unmatched_reads 7236844
Percentage reads assigned to barcode 81.5529
```
#### Matched sequences
Each assigned feature sequence is listed in the "feature_sequences.txt" file.
```

Feature Index Sequence Hamming_distance Counts Feature_name
  1 CAACTGCGTCCATGAAACAATAGACGCAGTTGAGAGTGGC  0      80 5_PDX1
  1 CAACTGCGTCCATGAAACAATAGACGCAGTTGAGAGaGGC  1       1 5_PDX1
  1 CAACTGaGTCCATGAAACAATAGACGCAGTTGAGAGTGGC  1       1 5_PDX1
  1 CAACTGCGTCCATaAAACAATAtACGCAGTTGAGAGTGGC  2       1 5_PDX1
  1 CAACTGCGTCCAgGAAACAATAGAaGCAGTTGAGAGTGGC  2       1 5_PDX1
  1 CAACTGatTaCATGAAACAATAGACGaAGTTtAGAGTGGC  5       1 5_PDX1
  2 GGTATGTGAACATACAACATAGGAGTTGGTTACAAGGAAT  0     109 12_PAX6
  2 GGTATGTGAACATACAACATAGGAGTTaGTTACAAGGAAT  1       2 12_PAX6
  2 GGTATGTGAACATACAcCATAGGAGTTGGTTACAAGGAAT  1       1 12_PAX6
  2 GGTATGTGAACATACAAaATAGGAGTTGGTTACAAGGAAT  1       1 12_PAX6
  2 GGTATGTGAACATACAACATAGGAtTaaGTTACAAGGAAT  3       1 12_PAX6
  3 TTGCTGACGACAGCTTACAGGCGAAACGGTCTTAAGTAAA  0      28 12_PAX6-alt1
  3 TTGCTGACGACgGCTTACAGGCGAAACGGTCTTAAGTAAA  1       2 12_PAX6-alt1
  4 ACGGTGCGAACATAAAACTAAAGACATAGTACTTAGATGT  0     233 15_NKX2-2
  4 ACGGTGCGAACATAAAACTAAAGACcTAGTACTTAGATGT  1       2 15_NKX2-2

```
Each of the matched sequences is displayed under the feature index that they are matched to.
#### Feature frequency histogram heatmap
After deduplication, we can make a histogram of the number of features assigned to cells. The perfect scenario would be nearly all cells have 1 feature assigned to them with a small distribution representing multiplets and another small distribution from empty cells (ambient counts that are assigned a cell barcode). Underneath is a heatmap where each row is the same histogram, but only showing cells that have that feature in it. The values of the histograms are space separated in feature_histograms.txt.
