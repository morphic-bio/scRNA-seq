# Morphic scRNA-seq 

### Organization of counts files
The files are organized by experiments (A,B,C,D,E,F,G_1,G_2,H,I,J,L_1,L_2)
For each of the the experiments there are two files for the different aligners (except CellRanger) containing the counts and feature counts (larry barcodes) which are organized as shown for experiment A
```
├── A
│   ├── expression_feature_assignment_QC
│   │   ├── feature_histograms.txt
│   │   ├── feature_sequences.txt
│   │   ├── heatmap.png
│   │   └── stats.txt
│   ├── kallisto
│   │   ├── unfiltered_counts.h5ad
│   │   └── unfiltered_features.h5ad
│   ├── piscem.splicei
│   │   ├── unfiltered_counts.h5ad
│   │   └── unfiltered_features.h5ad
│   ├── star
│   │   ├── unfiltered_counts.h5ad
│   │   └── unfiltered_features.h5ad
│   └── targeted_features_assignment_QC
│       ├── feature_histograms.txt
│       ├── feature_sequences.txt
│       ├── heatmap.png
│       └── stats.txt
├── B
   ...

```
### Structure of counts files
unfiltered_counts.h5ad
```
├── X (Data Matrix): (2035634, 36601)
├── obs: (2035634, 12)
│   └── is_cell (star)
│   └── non_empty
│   └── doublet
│   └── doublet_score
│   └── feature-10
│   └── feature-4
│   └── n_genes
│   └── mt_counts
│   └── mt_pct
│   └── filter
│   └── singlet
│   └── singlet_filtered
├── var: (36601, 0)
├── layers
│   └── ambiguous
│   └── denoised
│   └── spliced
│   └── unspliced
```
unfiltered_features.h5ad
```
├── X (Data Matrix): (2035634, 105)
├── obs: (2035634, 12)
│   └── is_cell (star)
│   └── non_empty
│   └── doublet
│   └── doublet_score
│   └── feature-10
│   └── feature-4
│   └── n_genes
│   └── mt_counts
│   └── mt_pct
│   └── filter
│   └── singlet
│   └── singlet_filtered
├── var: (105, 0)
├── layers
│   └── expression
│   └── sc_counts
```

### Definition of X, obs, var, layers
#### X (Data Matrix): (<number_of_barcodes>, 36601)
This is the basic 2D matrix with gene counts or feature counts. The gene counts include the feature counts as they are expressed. These are at the columns that follow the normal gene columns. The number of rows depends on the aligner. The gene counts include intron counts (as per CellRanger >= 7.0) and is the sum of the ambiguous, spliced and unspliced layers. The feature counts are determined using Biodepot's software (github under preparation).
#### obs: (<number_of_barcodes>, 12)
##### is_cell
Cellranger and Star-solo run their own empty-cell detection. This has a value of 1 when the internal methodology identifies the barcode to be a non-empty cell.
##### non_empty
The EmptyDropsCellRanger [algorithm](https://github.com/MarioniLab/DropletUtils/pull/64/commits/5a0b6c2f91b4767b92e3fb0894a2c05c9eaa22d4) is run to detect empty cells. A value of 1 indicates that the cell is non-empty
##### doublet
[scDblFinder](https://github.com/plger/scDblFinder) was run to detect doublets. A value of 1 indicates that the barcode is a doublet
##### doublet_score
Confidence that the cell is a doublet. 
##### feature-10
Name of feature if at least 10 feature counts found and second best no more than 5. Otherwise empty (Nan)
##### feature-4
Name of feature if at least 4 feature counts found and second best no more than 2. Otherwise empty (Nan)
#### n_genes
Number of genes with non-zero counts
##### mt_counts
Total number of counts from mitochondrial genes
##### mt_pct
Percentage of counts that are from mitochondrial genes
##### filter
Has value of 1 if  200 >= n_genes <= 2500, and mt_pct < 5%
##### singlet
Has a value of 1 if the barcode is not detected as a doublet and it is not empty
##### singlet_filtered
Has a value of 1 if the barcode is a singlet and passes the filter
#### var: (36601, 0) (for gene counts file)
##### unspliced, spliced, ambiguous
These are the different count values (layers) for the barcode, gene count matrix. The unspliced, spliced and ambiguous counts are defined by this figure from the kallisto preprint  https://www.biorxiv.org/content/10.1101/2022.12.02.518832v2.full with M being mature/spliced, A, ambiguous and N, nascent/unspliced. The definitions are those defined by Velocyto for all the aligners except Kallisto. Note that the workaround is to assign Ambiguous counts to spliced counts, in which case all aligners have the same definition.

![[Pasted image 20240829001543.png]]
##### denoised
The denoised layer is produced by cellBender a deep-learning package from Broad that detects and removes ambient counts. Note that cellBender uses an initial empty-cell detection method and only calculates the denoised counts for non-empty cells. The rest of the layer is zero-filled. Note that it is possible but unlikely that some non-empty cells detected by the EmptyDrops algorithm may produce zero counts when denoised. This is unlikely as cellBender intentionally uses a very generous initial criterion for filtering out empty cells.
#### For feature counts file
####  X (Data Matrix): (<number_of_barcodes, 105)
This is the basic 2D matrix with feature counts. The counts are determined using the Biodepot software from the targeted sequencing. UMI deduping proceeds by finding the feature that has the most reads assigned to a barcode-umi. Matches of up to 3 ambiguous bases and a max Hamming distance of 5 were considered as long as the match was unambiguous. Sequence barcodes are fuzzy matched using the CellRanger methodology. 
##### expression
The feature barcodes (larry barcodes) are expressed and the scRNA-seq fastqs can be searched for matches directly. The same matching methods were used except that de-duping was done using the standard RNA-seq strategy - all barcode-umi counts that mapped to a feature are reduced to 1.
##### sc_layer
This uses the alignment to the GFP-feature sequence to determine the counts. Previously, the feature sequence alone was used which gave reasonable counts only with Kallisto. The more complete sequence should work with star and piscem though this has not been tested.

### QC files

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
