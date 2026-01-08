# Load required library
library(DESeq2)

# ====================
# 1. SET YOUR INPUT FILES HERE
# ====================

# Path to your count matrix (CSV file with genes as rows, samples as columns)
count_matrix_path <- "DEGs/pseudobulk/counts_and_metadata/pseudobulk_counts_T-cells-pdac_vs_normal.csv"

# Path to your metadata file (CSV with 'condition' column and sample names as row names)
metadata_path <- "DEGs/pseudobulk/counts_and_metadata/pseudobulk_metadata_T-cells-pdac_vs_normal.csv"

# Output directory for results
output_dir <- "DEGs/pseudobulk"

# Define your comparison groups
group1 <- "PDAC"  # Experimental group
group2 <- "Normal"    # Control group

# DESeq2 parameters
lfc_threshold <- 1      # Log fold change threshold
alpha_value <- 0.05     # Adjusted p-value threshold
test_type <- "Wald"     # Statistical test: "Wald" or "LRT"
fit_type <- "parametric" # Fitting type: "parametric", "local", or "mean"
alt_hypothesis <- "greaterAbs" # Alternative hypothesis

# ====================
# 2. LOAD YOUR DATA
# ====================

# Read count matrix
count_matrix <- read.csv(count_matrix_path, row.names = 1, check.names = FALSE)
cat("Count matrix loaded:", nrow(count_matrix), "genes,", ncol(count_matrix), "samples\n")

# Read metadata
metadata <- read.csv(metadata_path, row.names = 1)
cat("Metadata loaded:", nrow(metadata), "samples\n")

# Check if 'condition' column exists
if (!"condition" %in% colnames(metadata)) {
  stop("ERROR: Metadata must have a 'condition' column")
}

# ====================
# 3. PREPARE DATA FOR DESeq2
# ====================

# Ensure samples match between count matrix and metadata
common_samples <- intersect(colnames(count_matrix), rownames(metadata))
cat("Common samples found:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("ERROR: No common samples between count matrix and metadata")
}

# Subset data to common samples (If you want to proceed with limited data)
#count_matrix <- count_matrix[, common_samples, drop = FALSE]
#metadata <- metadata[common_samples, , drop = FALSE]

# Ensure condition is a factor
metadata$condition <- factor(metadata$condition)

# ====================
# 4. RUN DESeq2 ANALYSIS
# ====================

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

# Keep genes with at least, e.g., 10 total counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2 with specified parameters
dds <- DESeq(
  dds,
  test = test_type,
  fitType = fit_type
)

# Get results for specific contrast
results <- results(
  dds,
  contrast = c("condition", group1, group2),
  alpha = alpha_value,
  lfcThreshold = lfc_threshold,
  altHypothesis = alt_hypothesis
)

# ====================
# 5. SAVE RESULTS
# ====================

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save full results
results_df <- as.data.frame(results)
write.csv(results_df, file.path(output_dir, "T-cells_pdac_vs_normal.csv"))
cat("Results saved to:", file.path(output_dir, "deseq2_results.csv"), "\n")

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file.path(output_dir, "normalized_counts.csv"))
cat("Normalized counts saved to:", file.path(output_dir, "normalized_counts.csv"), "\n")

# Save DESeq2 object for later use
saveRDS(dds, file.path(output_dir, "deseq2_object.rds"))
cat("DESeq2 object saved to:", file.path(output_dir, "deseq2_object.rds"), "\n")

# ====================
# 6. VIEW SUMMARY
# ====================

cat("\n========== ANALYSIS SUMMARY ==========\n")
cat("Comparison:", group1, "vs", group2, "\n")
cat("Total genes:", nrow(results_df), "\n")
cat("Significant genes (padj <", alpha_value, "):", sum(results_df$padj < alpha_value, na.rm = TRUE), "\n")
cat("Up-regulated (padj <", alpha_value, " & LFC > 0):", sum(results_df$padj < alpha_value & results_df$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("Down-regulated (padj <", alpha_value, " & LFC < 0):", sum(results_df$padj < alpha_value & results_df$log2FoldChange < 0, na.rm = TRUE), "\n")
cat("=====================================\n")

# View top 10 significant genes
cat("\nTop 10 significant genes:\n")
top_genes <- results_df[order(results_df$padj), ][1:10, ]
print(top_genes)
