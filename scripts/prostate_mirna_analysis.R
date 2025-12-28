# Prostate miRNA Differential Expression Analysis
# Template script for analyzing miRNA expression data from GEO dataset GSE60117
# This script performs a basic workflow to identify miRNAs differentially expressed
# between normal prostate tissue and prostate adenocarcinoma.
#
# Dataset: GSE60117 (miRNA expression profiling of prostate cancer and normal prostatic samples)
# Platform: Agilent-021827 Human miRNA Microarray (V3) (GPL13264)
# Samples: 21 normal, 56 tumor
#
# Author: Prostate-miRNA-Analysis-Starter
# Date: 2025-12-28

# ==================================================
# 0. Install required packages (if not already installed)
# ==================================================
# Uncomment the following lines to install packages:
# install.packages(c("BiocManager", "ggplot2", "dplyr", "tidyverse"))
# BiocManager::install(c("GEOquery", "limma"))

# Load libraries
library(GEOquery)
library(limma)
library(ggplot2)
library(dplyr)

# ==================================================
# 1. Download and load the dataset
# ==================================================
# Option 1: Download the series matrix file directly from GEO
# This may take a few minutes depending on your internet connection.
geo_accession <- "GSE60117"
gse <- getGEO(geo_accession, GSEMatrix = TRUE)

# The object 'gse' is a list; the expression matrix is usually in the first element
if (length(gse) > 0) {
  expr_data <- exprs(gse[[1]])
  pheno_data <- pData(phenoData(gse[[1]]))
} else {
  stop("Failed to retrieve data from GEO.")
}

# Option 2: If you have already downloaded the series matrix file locally,
# load it using:
# gse <- getGEO(filename = "GSE60117_series_matrix.txt.gz")

# Inspect dimensions
cat("Expression matrix dimensions:", dim(expr_data), "\n")
cat("Phenotype data dimensions:", dim(pheno_data), "\n")

# ==================================================
# 2. Preprocessing and normalization
# ==================================================
# The data from GEO may already be log2 transformed and normalized.
# Check a few values to decide whether further normalization is needed.
summary(as.vector(expr_data[1:10, 1:5]))

# If the data are not logged, apply log2 transformation (skip if already logged)
# expr_data <- log2(expr_data + 1)

# If needed, perform between‑array normalization (e.g., quantile normalization)
# expr_data_normalized <- normalizeBetweenArrays(expr_data, method = "quantile")
# expr_data <- expr_data_normalized

# ==================================================
# 3. Prepare phenotype groups
# ==================================================
# The phenotype table contains sample characteristics.
# Identify which samples are normal and which are tumor.
# In GSE60117, the column "source_name_ch1" or "characteristics_ch1" may contain this information.
head(pheno_data[, 1:5])

# For this dataset, we can use the column "source_name_ch1" to define groups.
# Adjust the column name if different.
if ("source_name_ch1" %in% colnames(pheno_data)) {
  pheno_data$group <- ifelse(grepl("normal", pheno_data$source_name_ch1, ignore.case = TRUE),
                             "Normal", "Tumor")
} else {
  # If the column is not present, you may need to inspect other columns or
  # manually assign groups based on sample names.
  stop("Please inspect pheno_data to define groups. Modify the script accordingly.")
}

# Create a factor for the design matrix
group <- factor(pheno_data$group, levels = c("Normal", "Tumor"))
table(group)

# ==================================================
# 4. Differential expression analysis with limma
# ==================================================
# Create design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- c("Normal", "Tumor")

# Fit linear model
fit <- lmFit(expr_data, design)

# Define contrast (Tumor vs Normal)
contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract results
results <- topTable(fit2, number = Inf, adjust.method = "BH")
head(results)

# ==================================================
# 5. Save results
# ==================================================
# Create results directory if it doesn't exist
if (!dir.exists("../results")) {
  dir.create("../results", recursive = TRUE)
}

# Write full results to a CSV file
write.csv(results, file = "../results/differential_expression_results.csv", row.names = TRUE)

# Filter for significant miRNAs (adjusted p-value < 0.05 and |logFC| > 1)
significant <- results %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

cat("Number of significant miRNAs:", nrow(significant), "\n")
write.csv(significant, file = "../results/significant_miRNAs.csv", row.names = TRUE)

# ==================================================
# 6. Visualization: Volcano plot
# ==================================================
# Create a volcano plot using ggplot2
volcano_data <- results %>%
  mutate(Significance = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                               ifelse(logFC > 0, "Upregulated", "Downregulated"),
                               "Not significant"))

# Create the plot
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated" = "blue",
                                "Upregulated" = "red",
                                "Not significant" = "gray60")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(title = "Volcano plot: miRNA expression in prostate cancer vs normal",
       x = "log2 Fold Change (Tumor vs Normal)",
       y = "-log10(adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the plot
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}
ggsave(filename = "../figures/volcano_plot.png", plot = volcano_plot,
       width = 8, height = 6, dpi = 300)

# Display the plot
print(volcano_plot)

# ==================================================
# 7. Session info
# ==================================================
# Record R session information for reproducibility
sessionInfo()