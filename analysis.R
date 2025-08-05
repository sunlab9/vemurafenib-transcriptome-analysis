# Load necessary libraries
library(GEOquery)                          # To access data from the GEO database
library(affy)                              # For reading and processing Affymetrix microarray data
library(hugene10stv1cdf)                   # CDF environment for the Affymetrix Human Gene 1.0 ST Array
library(hugene10sttranscriptcluster.db)    # Annotation database for the chip
library(AnnotationDbi)                     # Database interface for annotation
library(org.Hs.eg.db)                      # Gene annotations for Homo sapiens
library(limma)                             # Linear models for microarray data
library(pheatmap)                          # For creating heatmaps
library(ggplot2)                           # For data visualization
library(clusterProfiler)                   # Functional enrichment analysis
library(EnhancedVolcano)                   # For enhanced volcano plots
library(pathview)                          # For KEGG pathway visualization

# Unzip the raw microarray data
untar("GSE42872_RAW.tar", exdir = "GSE42872_RAW")

# Read CEL files (raw Affymetrix data)
data <- ReadAffy(celfile.path = "GSE42872_RAW")

# Normalize the data using RMA (Robust Multi-array Average)
norm_data <- rma(data)

# Extract expression matrix
exprs_matrix <- exprs(norm_data)

# Map probe IDs to gene symbols
gene_symbols <- mapIds(hugene10sttranscriptcluster.db,
                       keys = rownames(exprs_matrix),
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Combine expression data with gene symbols
exprs_annotated <- data.frame(GeneSymbol = gene_symbols, exprs_matrix)

# Remove rows with missing gene symbols
exprs_clean <- exprs_annotated[!is.na(exprs_annotated$GeneSymbol), ]

# Remove duplicated gene symbols, keeping the first occurrence
exprs_dedup <- exprs_clean[!duplicated(exprs_clean$GeneSymbol), ]

# Set gene symbols as row names and remove GeneSymbol column
rownames(exprs_dedup) <- exprs_dedup$GeneSymbol
exprs_dedup$GeneSymbol <- NULL

# Rename columns for clarity
colnames(exprs_dedup) <- c("Control_rep1", "Control_rep2", "Control_rep3",
                           "Vemurafenib_rep1", "Vemurafenib_rep2", "Vemurafenib_rep3")

# Define experimental groups
group <- factor(c(rep("Control", 3), rep("Vemurafenib", 3)))

# Create design matrix for differential expression analysis
design <- model.matrix(~ group)

# Fit linear model to the expression data
fit <- lmFit(exprs_dedup, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract all differentially expressed genes
deg_results <- topTable(fit, coef = "groupVemurafenib", number = Inf, adjust = "fdr")

# Extract log fold changes and adjusted p-values
logFC <- deg_results$logFC
adj.P.val <- deg_results$adj.P.Val

# Create a volcano plot to visualize DEGs
plot(logFC, -log10(adj.P.val),
     pch = 20,
     main = "Volcano Plot",
     xlab = "log2 Fold Change",
     ylab = "-log10 Adjusted P-value",
     col = ifelse(adj.P.val < 0.05 & abs(logFC) > 1, "red", "black"))

# Add cutoff lines to the volcano plot
abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)

# Get top 10 significant DEGs
top10 <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 1)[1:10, ]

# Save top 10 DEGs to a CSV file
write.csv(top10, file = "top_10_deg_genes.csv", row.names = TRUE)

# Extract expression matrix for top 10 DEGs
expr_top10 <- exprs_dedup[rownames(top10), ]

# Scale data for heatmap
expr_top10_scaled <- t(scale(t(expr_top10)))

# Create annotation for samples
annotation_df <- data.frame(Group = group)
rownames(annotation_df) <- colnames(expr_top10_scaled)

# Plot heatmap of top 10 DEGs
pheatmap(expr_top10_scaled,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         fontsize = 10)

# Filter significant DEGs for enrichment analysis
deg_filtered <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 1)
deg_filtered <- deg_filtered[order(deg_filtered$adj.P.Val), ]

# Select top 100 genes for GO/KEGG
top_genes <- head(rownames(deg_filtered), 100)

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(top_genes,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# Perform GO enrichment for Biological Process (BP)
ego_BP <- enrichGO(gene = gene_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

# Plot GO BP results
barplot(ego_BP, showCategory = 10, title = "GO Biological Process")

# Perform GO enrichment for Molecular Function (MF)
ego_MF <- enrichGO(gene = gene_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

# Plot GO MF results
barplot(ego_MF, showCategory = 10, title = "GO Molecular Function")

# Perform GO enrichment for Cellular Component (CC)
ego_CC <- enrichGO(gene = gene_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

# Plot GO CC results
barplot(ego_CC, showCategory = 10, title = "GO Cellular Component")

# Prepare gene fold change vector for KEGG analysis
gene_fc <- deg_results$logFC
names(gene_fc) <- rownames(deg_results)

# Map gene symbols to Entrez IDs for KEGG
gene_fc_entrez <- bitr(names(gene_fc),
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)

# Align fold changes with Entrez IDs
gene_fc_mapped <- gene_fc[gene_fc_entrez$SYMBOL]
names(gene_fc_mapped) <- gene_fc_entrez$ENTREZID

# Visualize gene expression on KEGG Melanoma pathway
pathview(gene.data  = gene_fc_mapped,
         pathway.id = "hsa05200",  # Melanoma pathway ID
         species    = "hsa")








