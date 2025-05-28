# Metabolic genes expression of E. coli under the experimental evolution assay with different antibiotics

## Project Description
This project analyzes transcriptomic data from experimental evolution studies of E. coli under antibiotic selection pressure. The workflow processes data from multiple studies to identify consistent gene expression changes associated with antibiotic resistance development.

## Pipeline Components

1. **Data Preprocessing**
   - Merges data from different studies
   - Filters for metabolic genes of interest
   - Handles missing values and normalization

2. **Differential Expression Analysis**
   - Computes log2 fold changes between treated and control samples
   - Calculates p-values and adjusted q-values
   - Combines results across studies using weighted averages

3. **Visualization**
   - Volcano plots showing significant expression changes
   - Heatmaps of gene expression patterns across antibiotics

## Requirements

### Python
- pandas
- numpy

### R
- tidyverse (dplyr, tidyr, ggplot2)
- ComplexHeatmap
- circlize
- ggrepel

## Usage

1. Place your input files in `data/raw/`:
   - `metabolic_genes.csv` - Gene list of interest
   - `file.csv` - Transcriptomic data
   - Additional study files as needed

2. Run the pipeline in order:
   - `1_data_filtering.py` (Python)
   - `2_analysis.R` (R)
   - `3_visualization.R` (R)

3. Outputs will be saved in `results/`:
   - Filtered data files
   - Statistical results
   - Publication-quality figures

# File Structure


├── data/
│   ├── data1.csv
│   ├── data2.csv
│   └── ...
├── scripts/
│   ├── filter_genes.py
│   ├── merge_data.R
│   └── visualization.R
├── output/
│   ├── Filtered data.csv
│   ├── combined_data_wide_format.csv
│   ├── volcano_plot_high_res.png
│   ├── heatmap_significant_genes_high_res.png
│   └── gene_expression_per_antibiotic_columns.csv
├── README.md
└── ...
Usage Instructions
1. Filtering Genes (Python)


python scripts/filter_genes.py
Input required:
data1.csv (full transcriptomic profiles)
data2.csv (metabolic gene list)
Output:
Filtered gene file: Filtered data.csv
2. Merging and Expression Analysis (R)

Rscript scripts/merge_data.R
Input required:
Filtered data files from multiple studies (e.g., filtered_genes1.csv, filtered_genes2.csv)
Output:
Merged data: combined_data_wide_format.csv
Per-antibiotic summary: gene_expression_per_antibiotic_columns.csv
3. Visualization (R)

Rscript scripts/visualization.R
Output:
Volcano and heatmap plots in output/ directory

## Requirements
Python 3.7+, with pandas
R 4.0+, with packages: dplyr, ggplot2, pheatmap, ComplexHeatmap, tidyr, clusterProfiler, ggrepel, circlize

## Install R dependencies:

install.packages(c("dplyr", "ggplot2", "pheatmap", "tidyr", "ggrepel"))
BiocManager::install(c("ComplexHeatmap", "clusterProfiler", "circlize"))
Important Notes
The provided pipeline assumes certain column names in the CSV files. Check your headers and adjust scripts accordingly!
Statistical calculations (p- and q-values) are based on a fixed variance. For robust differential expression analysis, use replicates and modeling tools like DESeq2 or edgeR.
Plots highlight significant genes with |log2FC| above threshold and p < 0.05.






