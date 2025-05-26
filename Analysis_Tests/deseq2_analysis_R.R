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


# Create ranked gene list
# Sort the DEGs by log2FoldChange (descending)
ranked_genes <- res[order(res$log2FoldChange, decreasing = TRUE), ]

# Create a data frame with gene names and log2FoldChange values
ranked_gene_list <- data.frame(
  gene = rownames(ranked_genes),
  log2FoldChange = ranked_genes$log2FoldChange
)

# Write the ranked gene list to a file for GSEA
write.table(ranked_gene_list, file = "ranked_gene_list.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# View the first few rows of the ranked gene list
head(ranked_gene_list)
