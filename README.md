# Prostate-miRNA-Analysis-Starter

## Overview
This repository provides a practical, educational guide for performing a basic differential expression analysis of miRNAs to distinguish between normal and malignant prostate tissue. It serves as a "starter kit" for researchers or students entering the field of cancer biomarker discovery using miRNA expression data.

The analysis workflow covers data loading, preprocessing, normalization, differential expression testing, and visualization using a publicly available dataset from the NCBI Gene Expression Omnibus (GEO).

## Dataset
The analysis uses the GEO dataset **GSE60117**:
- **Title**: miRNA expression profiling of prostate cancer (PCa) and normal prostatic samples
- **Platform**: Agilent-021827 Human miRNA Microarray (V3) (GPL13264)
- **Samples**: 21 normal prostate tissues, 56 prostate cancer tissues
- **GEO Accession**: [GSE60117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60117)
- **Publication**: PMID: 27384993

## Repository Structure
```
Prostate-miRNA-Analysis-Starter/
├── README.md               # This file
├── scripts/
│   └── prostate_mirna_analysis.R   # R script for differential expression analysis
├── results/                # Directory for output tables
└── figures/                # Directory for generated plots
```

## Getting Started

### Prerequisites
- **R** (version 4.0 or later) or **RStudio**
- Required R packages: `limma`, `ggplot2`, `dplyr`, `tidyverse` (optional)

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/britdev071809/Prostate-miRNA-Analysis-Starter.git
   cd Prostate-miRNA-Analysis-Starter
   ```
2. Install the required R packages:
   ```r
   install.packages(c("limma", "ggplot2", "dplyr"))
   ```

### Running the Analysis
1. Open the script `scripts/prostate_mirna_analysis.R` in R/RStudio.
2. Follow the step‑by‑step comments to download the dataset, preprocess the data, perform differential expression analysis, and generate a volcano plot.
3. Output tables will be saved in the `results/` directory, and plots in the `figures/` directory.

## Analysis Steps
The script outlines the following workflow:
1. **Data acquisition**: Download the GEO dataset using GEOquery or load the provided expression matrix.
2. **Data preprocessing**: Log‑transform (if needed), normalize between arrays (e.g., quantile normalization), and filter low‑expressed miRNAs.
3. **Differential expression**: Use the `limma` package to fit a linear model and compute moderated t‑statistics, p‑values, and log‑fold‑changes between normal and cancer samples.
4. **Visualization**: Generate a volcano plot highlighting significantly up‑ and down‑regulated miRNAs.
5. **Output**: Save a table of differentially expressed miRNAs (e.g., with adjusted p‑value < 0.05 and |logFC| > 1).

## Expected Results
The analysis will identify a set of miRNAs that are differentially expressed between normal prostate tissue and prostate adenocarcinoma. These candidate biomarkers can serve as a starting point for further validation or functional studies.

## Contributing
Contributions are welcome! If you have suggestions for improving the analysis, adding alternative scripts (e.g., in Python), or extending the documentation, please open an issue or submit a pull request.

## License
This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

## Citation
If you use this tutorial or the dataset in your work, please cite the original publication:
- **GSE60117**: PMID 27384993
a