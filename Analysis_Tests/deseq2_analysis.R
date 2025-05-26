# Load DESeq2
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
library(DESeq2)

# Load ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")
}
library(ggplot2)


# Define file paths
files <- c(
  "600" = "/scratch/general/vast/u0962361/working/counts_STAR/600.tab",
  "601" = "/scratch/general/vast/u0962361/working/counts_STAR/601.tab",
  "602" = "/scratch/general/vast/u0962361/working/counts_STAR/602.tab"
)

# Read STAR ReadsPerGene.out.tab (column 2 = unstranded counts)
read_counts <- lapply(files, function(path) {
  dat <- read.table(path, header = FALSE)
  dat <- dat[, c(1, 2)]  # Column 1 = gene name, Column 2 = unstranded count
  colnames(dat) <- c("gene", "count")
  return(dat)
})

# Merge count data
merged <- Reduce(function(x, y) merge(x, y, by = "gene"), read_counts)
rownames(merged) <- merged$gene
merged <- merged[!grepl("^__", merged$gene), ]  # Remove STAR summary rows
count_matrix <- merged[, -1]
colnames(count_matrix) <- names(files)

# Sample metadata
col_data <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(c("control", "control", "treatment"))
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~condition)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)

# Save results
write.csv(as.data.frame(res), "deseq2_results.csv")

# Filter significant DEGs
sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# Order by adjusted p-value (most significant first)
top_genes <- sig_genes[order(sig_genes$padj), ]

# Get top 10 gene names
top_labels <- rownames(head(top_genes, 10))


# MA Plot
pdf("MA_plot.pdf")
plotMA(res, ylim=c(-5,5))
dev.off()

# Convert results to data frame and add gene names
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Load ggrepel if not already installed
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)

# Volcano plot with top gene labels
pdf("volcano_plot_labeled.pdf")
ggplot(res_df, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.5) +
  geom_text_repel(data = subset(res_df, gene %in% top_labels),
                  aes(label = gene), size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(p-value)")
dev.off()

# Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)

# PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Plot PCA
pdf("PCA_plot.pdf")
ggplot(pca_data, aes(PC1, PC2, color = condition)) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(100 * pca_data$percentVar[1]), "% variance")) +
  ylab(paste0("PC2: ", round(100 * pca_data$percentVar[2]), "% variance")) +
  theme_minimal() +
  ggtitle("PCA Plot")

dev.off()
