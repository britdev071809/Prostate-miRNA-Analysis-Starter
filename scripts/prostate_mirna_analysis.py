#!/usr/bin/env python3
"""
Prostate miRNA Differential Expression Analysis (Python version)
Template script for analyzing miRNA expression data from GEO dataset GSE60117.
This script performs a basic workflow to identify miRNAs differentially expressed
between normal prostate tissue and prostate adenocarcinoma.

Dataset: GSE60117 (miRNA expression profiling of prostate cancer and normal prostatic samples)
Platform: Agilent-021827 Human miRNA Microarray (V3) (GPL13264)
Samples: 21 normal, 56 tumor

Author: Prostate-miRNA-Analysis-Starter
Date: 2025-12-28
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ==================================================
# 0. Install required packages (if not already installed)
# ==================================================
# Run the following commands in your terminal:
# pip install pandas numpy scipy statsmodels matplotlib seaborn
# Optional: pip install GEOparse (for direct GEO download)

# ==================================================
# 1. Load the dataset
# ==================================================
# Option 1: Use GEOparse to download and parse the dataset directly
# Uncomment the following lines if you have GEOparse installed:
"""
import GEOparse
gse = GEOparse.get_GEO(geo="GSE60117", destdir="./")
expr_data = gse.pivot_samples('VALUE')  # expression matrix (genes x samples)
pheno_data = gse.phenotype_data         # sample metadata
"""

# Option 2: Load a pre‑downloaded series matrix file (recommended for beginners)
# Download the series matrix file from:
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60117&format=file
# Save it as 'GSE60117_series_matrix.txt.gz' in the same directory as this script.
# The following code reads the series matrix using pandas.

def load_series_matrix(filepath):
    """Parse a GEO series matrix file and return expression and phenotype data."""
    # This is a simplified parser; in practice you may need to adjust based on file format.
    # For a robust solution consider using GEOparse.
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find start of expression matrix (usually after "!series_matrix_table_begin")
    start_idx = None
    for i, line in enumerate(lines):
        if line.startswith('!series_matrix_table_begin'):
            start_idx = i + 1
            break
    if start_idx is None:
        raise ValueError("Could not find expression table in series matrix file.")
    
    # The line after the begin marker contains column headers (sample IDs)
    header = lines[start_idx].strip().split('\t')
    # The first column is the probe/gene identifier
    # Read the data rows until "!series_matrix_table_end"
    data_rows = []
    for i in range(start_idx + 1, len(lines)):
        if lines[i].startswith('!series_matrix_table_end'):
            break
        data_rows.append(lines[i].strip().split('\t'))
    
    expr_df = pd.DataFrame(data_rows, columns=header)
    expr_df.set_index(header[0], inplace=True)
    # Convert to numeric
    expr_df = expr_df.apply(pd.to_numeric, errors='coerce')
    
    # Extract phenotype information from the header lines before the table
    # This is a minimal extraction; you may need to parse more carefully.
    pheno_lines = [line for line in lines if line.startswith('!Sample_')]
    pheno_dict = {}
    for line in pheno_lines:
        parts = line.strip().split(' = ')
        if len(parts) == 2:
            key = parts[0].replace('!Sample_', '')
            val = parts[1]
            pheno_dict.setdefault(key, []).append(val)
    pheno_df = pd.DataFrame(pheno_dict)
    
    return expr_df, pheno_df

# Uncomment and adjust the file path as needed
# expr_data, pheno_data = load_series_matrix('GSE60117_series_matrix.txt.gz')

# For the purpose of this template, we will simulate a dummy expression matrix
# and phenotype labels. Replace this with actual data loading.
print("WARNING: Using simulated data for demonstration. Replace with actual data loading.")
# Simulate expression data (100 miRNAs, 77 samples)
np.random.seed(42)
expr_data = pd.DataFrame(np.random.randn(100, 77),
                         index=[f'miR-{i}' for i in range(1, 101)],
                         columns=[f'Sample_{i}' for i in range(1, 78)])
# Simulate groups: first 21 samples normal, rest tumor
groups = ['Normal'] * 21 + ['Tumor'] * 56
pheno_data = pd.DataFrame({'sample': expr_data.columns, 'group': groups})

# ==================================================
# 2. Preprocessing and normalization
# ==================================================
# The data may already be log2 transformed; check distribution.
# If not logged, apply log2 transformation (skip if already logged).
# expr_data = np.log2(expr_data + 1)

# Optional: quantile normalization (requires sklearn or other packages)
# from sklearn.preprocessing import quantile_transform
# expr_data_normalized = pd.DataFrame(
#     quantile_transform(expr_data.T, n_quantiles=min(expr_data.shape[1], 1000),
#                        output_distribution='normal').T,
#     index=expr_data.index, columns=expr_data.columns)
# expr_data = expr_data_normalized

# ==================================================
# 3. Prepare phenotype groups
# ==================================================
# In real data, use pheno_data to define groups.
# For the simulated data we already have groups.
group_series = pheno_data.set_index('sample')['group']
# Ensure sample order matches expression matrix
group_series = group_series.loc[expr_data.columns]
print("Group counts:")
print(group_series.value_counts())

# ==================================================
# 4. Differential expression analysis (t-test)
# ==================================================
# Separate expression matrices for normal and tumor
normal_expr = expr_data.loc[:, group_series == 'Normal']
tumor_expr = expr_data.loc[:, group_series == 'Tumor']

# Perform independent t-test for each miRNA
results = []
for mirna in expr_data.index:
    normal_vals = normal_expr.loc[mirna].values
    tumor_vals = tumor_expr.loc[mirna].values
    t_stat, p_val = stats.ttest_ind(tumor_vals, normal_vals, equal_var=False)  # Welch's t-test
    logfc = np.mean(tumor_vals) - np.mean(normal_vals)  # simple mean difference
    results.append({
        'miRNA': mirna,
        'logFC': logfc,
        't_stat': t_stat,
        'p_value': p_val
    })

results_df = pd.DataFrame(results)
results_df.set_index('miRNA', inplace=True)

# Adjust p-values for multiple testing (Benjamini-Hochberg FDR)
results_df['adj_p_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]

# Sort by adjusted p-value
results_df.sort_values('adj_p_value', inplace=True)
print("Top differentially expressed miRNAs:")
print(results_df.head(10))

# ==================================================
# 5. Save results
# ==================================================
# Create results directory if it doesn't exist
os.makedirs('../results', exist_ok=True)

# Write full results to a CSV file
results_df.to_csv('../results/differential_expression_results_python.csv')

# Filter for significant miRNAs (adj. p-value < 0.05 and |logFC| > 1)
significant = results_df[(results_df['adj_p_value'] < 0.05) & (abs(results_df['logFC']) > 1)]
print(f"Number of significant miRNAs: {significant.shape[0]}")
significant.to_csv('../results/significant_miRNAs_python.csv')

# ==================================================
# 6. Visualization: Volcano plot
# ==================================================
# Prepare data for volcano plot
volcano_data = results_df.copy()
volcano_data['-log10(p)'] = -np.log10(volcano_data['adj_p_value'])
volcano_data['Significance'] = 'Not significant'
volcano_data.loc[(volcano_data['adj_p_value'] < 0.05) & (volcano_data['logFC'] > 1), 'Significance'] = 'Upregulated'
volcano_data.loc[(volcano_data['adj_p_value'] < 0.05) & (volcano_data['logFC'] < -1), 'Significance'] = 'Downregulated'

# Plot
plt.figure(figsize=(8, 6))
sns.scatterplot(data=volcano_data, x='logFC', y='-log10(p)', hue='Significance',
                palette={'Downregulated': 'blue', 'Upregulated': 'red', 'Not significant': 'gray'},
                alpha=0.6, s=40)
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)
plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
plt.title('Volcano plot: miRNA expression in prostate cancer vs normal (Python)')
plt.xlabel('log2 Fold Change (Tumor vs Normal)')
plt.ylabel('-log10(adjusted p-value)')
plt.legend(title='Significance', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/volcano_plot_python.png', dpi=300)
plt.show()

# ==================================================
# 7. Session info
# ==================================================
print("\nPython packages used:")
import sys
print(f"Python {sys.version}")
for pkg in [pd, np, stats, plt]:
    print(f"{pkg.__name__} {pkg.__version__}")