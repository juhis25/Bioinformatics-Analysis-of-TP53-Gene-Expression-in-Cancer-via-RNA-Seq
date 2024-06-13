# Load libraries
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(pheatmap)

# Load count data and metadata
countData <- read.csv("processed_gene_expression_data.csv")
metadata <- read.csv("metadata.csv")

# Prepare DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countData[, 4:7]), colData = metadata, design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Filter significant genes
sig_genes <- res[which(res$padj < 0.05), ]

# Volcano plot
volcano_plot <- ggplot(data=as.data.frame(res), aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=padj < 0.05)) +
  theme_minimal() +
  labs(title="Volcano Plot of Differentially Expressed Genes", x="Log2 Fold Change", y="-Log10 Adjusted P-value")
print(volcano_plot)

# Heatmap of top 20 differentially expressed genes
top20_genes <- head(order(res$padj), 20)
heatmap_data <- assay(rlog(dds))[top20_genes, ]
pheatmap(heatmap_data, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, show_colnames=TRUE)


