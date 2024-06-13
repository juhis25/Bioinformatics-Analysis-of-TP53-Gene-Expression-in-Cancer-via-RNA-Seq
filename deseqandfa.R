# Prepare count matrix and metadata
countData <- as.matrix(data[, c("unstranded", "stranded_first", "stranded_second")])
colData <- metadata

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# View results
head(res)

# Normalized counts
normCounts <- counts(dds, normalized=TRUE)

# Plot a heatmap of the top differentially expressed genes
topGenes <- head(order(res$padj), 20)
heatmap(normCounts[topGenes, ])

# Volcano plot
ggplot(data=as.data.frame(res), aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=padj < 0.05)) +
  theme_minimal() +
  labs(title="Volcano Plot of Differentially Expressed Genes",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value")

# Prepare gene list
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = TRUE)

# Run GSEA
gseaResults <- gseKEGG(geneList = geneList, organism = 'hsa', pvalueCutoff = 0.05)
head(gseaResults)
