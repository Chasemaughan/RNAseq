library(DESeq2)

# Read counts
counts <- read.delim("data/feature_counts/ReadsPerGene.txt", comment.char = "#")

# Extract gene IDs and counts only
count_data <- counts[, c("Geneid", grep("\\.bam$", colnames(counts), value = TRUE))]
rownames(count_data) <- count_data$Geneid
count_data <- count_data[, -1]

# Rename columns (optional, but keeps things clean)
colnames(count_data) <- c("SRR24289600", "SRR24289601", "SRR24289602")

# Sample metadata
condition <- factor(c("control", "control", "treated"))
coldata <- data.frame(row.names = colnames(count_data), condition)

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Output results to file
write.csv(as.data.frame(res), file = "deseq2_results.csv")