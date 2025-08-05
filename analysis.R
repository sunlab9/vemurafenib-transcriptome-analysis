library(GEOquery)
library(affy)
library(hugene10stv1cdf)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(limma)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(EnhancedVolcano)
library(pathview)

untar("GSE42872_RAW.tar", exdir = "GSE42872_RAW")
data <- ReadAffy(celfile.path = "GSE42872_RAW")
norm_data <- rma(data)
exprs_matrix <- exprs(norm_data)


gene_symbols <- mapIds(hugene10sttranscriptcluster.db,
                       keys = rownames(exprs_matrix),
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")


exprs_annotated <- data.frame(GeneSymbol = gene_symbols, exprs_matrix)
exprs_clean <- exprs_annotated[!is.na(exprs_annotated$GeneSymbol), ]
exprs_dedup <- exprs_clean[!duplicated(exprs_clean$GeneSymbol), ]
rownames(exprs_dedup) <- exprs_dedup$GeneSymbol
exprs_dedup$GeneSymbol <- NULL

colnames(exprs_dedup) <- c("Control_rep1", "Control_rep2", "Control_rep3",
                           "Vemurafenib_rep1", "Vemurafenib_rep2", "Vemurafenib_rep3")




group <- factor(c(rep("Control", 3), rep("Vemurafenib", 3)))
design <- model.matrix(~ group)
fit <- lmFit(exprs_dedup, design)
fit <- eBayes(fit)
deg_results <- topTable(fit, coef = "groupVemurafenib", number = Inf, adjust = "fdr")

logFC <- deg_results$logFC
adj.P.val <- deg_results$adj.P.Val


plot(logFC, -log10(adj.P.val),
     pch = 20,
     main = "Volcano Plot",
     xlab = "log2 Fold Change",
     ylab = "-log10 Adjusted P-value",
     col = ifelse(adj.P.val < 0.05 & abs(logFC) > 1, "red", "black"))



abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)



top10 <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 1)[1:10, ]
write.csv(top10, file = "top_10_deg_genes.csv", row.names = TRUE)


expr_top10 <- exprs_dedup[rownames(top10), ]
expr_top10_scaled <- t(scale(t(expr_top10)))
annotation_df <- data.frame(Group = group)
rownames(annotation_df) <- colnames(expr_top10_scaled)

pheatmap(expr_top10_scaled,
         annotation_col = annotation_df,
         show_rownames = TRUE,
         fontsize = 10)

deg_filtered <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 1)
deg_filtered <- deg_filtered[order(deg_filtered$adj.P.Val), ]
top_genes <- head(rownames(deg_filtered), 100)

gene_entrez <- bitr(top_genes,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)


ego_BP <- enrichGO(gene = gene_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)
barplot(ego_BP, showCategory = 10, title = "GO Biological Process")


ego_MF <- enrichGO(gene = gene_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pvalueCutoff = 0.05,
                   readable = TRUE)
barplot(ego_MF, showCategory = 10, title = "GO Molecular Function")

ego_CC <- enrichGO(gene = gene_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pvalueCutoff = 0.05,
                   readable = TRUE)
barplot(ego_CC, showCategory = 10, title = "GO Cellular Component")

gene_fc <- deg_results$logFC
names(gene_fc) <- rownames(deg_results)

gene_fc_entrez <- bitr(names(gene_fc),
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)

gene_fc_mapped <- gene_fc[gene_fc_entrez$SYMBOL]
names(gene_fc_mapped) <- gene_fc_entrez$ENTREZID

pathview(gene.data  = gene_fc_mapped,
         pathway.id = "hsa05200",  # Melanoma pathway örneği
         species    = "hsa")








