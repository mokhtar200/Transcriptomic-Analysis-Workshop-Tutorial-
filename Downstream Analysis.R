# Step 1: Load Required Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "pheatmap", "ggplot2"))

library(DESeq2)
library(pheatmap)
library(ggplot2)

# Step 2: Import Kallisto Abundance Files
control_abundance <- read.delim("D:/Workshop/Control , 12Hr , Root/abundance.tsv", header = TRUE)
stress_abundance  <- read.delim("D:/Workshop/Stress, 12Hr, Root/abundance.tsv", header = TRUE)

# Combine estimated counts
tx_counts <- data.frame(
  Control = control_abundance$est_counts,
  Stress = stress_abundance$est_counts
)
rownames(tx_counts) <- control_abundance$target_id

tx_count_1=round(tx_counts)

# Step 3: Create Metadata
coldata <- data.frame(
  condition = factor(c("Control", "Stress"))
)

# Step 4: DESeq2 Analysis
# Adjust for lack of replicates
dds <- DESeqDataSetFromMatrix(countData = tx_count_1, 
                              colData = coldata, 
                              design = ~ condition)

# Estimate size factors and dispersion manually
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType = "local")

# Differential expression analysis
dds <- DESeq(dds, fitType = "local")
res <- results(dds, contrast = c("condition", "Stress", "Control"))

# Filter significant DEGs (adjusted p-value < 0.05)
res_sig <- subset(res, padj < 0.05)

# Save results
write.csv(as.data.frame(res), "DEG_results.csv")
write.csv(as.data.frame(res_sig), "Significant_DEGs.csv")

# Step 5: Visualization
# MA Plot
plotMA(res, main = "MA Plot (Stress vs Control)", ylim = c(-4, 4))

# Volcano Plot
res_df <- as.data.frame(res)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")

# Heatmap of Top 20 DEGs
vsd <- vst(dds, blind = FALSE)
top_genes <- head(order(res$padj), 20)
pheatmap(assay(vsd)[top_genes, ], cluster_rows = TRUE, show_rownames = TRUE,
         annotation_col = coldata)

# Step 6: Export Normalized Counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "Normalized_Counts.csv")
