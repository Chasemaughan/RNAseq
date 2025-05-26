library(DESeq2)         # For DESeqDataSet, DESeq, vst, results, plotMA
library(dplyr)          # For select (from tidyverse)
library(tibble)         # For column_to_rownames
library(pheatmap)       # For the heatmap
library(EnhancedVolcano) # For volcano plot

# Read counts
counts <- read.delim("data/feature_counts/ReadsPerGene.txt", comment.char = "#")

# Extract gene IDs and counts only
count_data <- counts %>%
  dplyr::select(Geneid, ends_with(".bam")) %>%
  column_to_rownames("Geneid")

# Rename columns for clarity
colnames(count_data) <- c("SRR24289600", "SRR24289601", "SRR24289602")

# Set up sample metadata
condition <- factor(c("control", "control", "treated"))
coldata <- data.frame(row.names = colnames(count_data), condition)


dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

top_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)
pheatmap(assay(vsd)[top_var_genes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE)


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = "Treated (SRR24289602) vs Control (SRR24289600/SRR24289601)",
                subtitle = "Differential Expression")


plotMA(res, ylim = c(-5, 5))
