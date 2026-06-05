library(DESeq2)
library(sva)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)

#Input_disease_vs_control_combined_count_file

count_data <- read.csv("TARGETP2vsControl_merged.csv", row.names = 1)

#remove-low_count_genes

count_data <- count_data[which(rowSums(count_data) > 100),]

#SegregateDisease_and_Control_samples

condition <- factor(c(rep("B-ALL", 198), rep("Control", 24)))
batch <- factor(c(rep("B1", 198), rep("B2", 24)))
coldata <- data.frame(row.names = colnames(count_data), condition, batch)

#Create_dds_object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ condition)

vsd <- vst(dds, blind=TRUE)
mat <- assay(vsd)

#PCA_before_SVA

pca_before <- prcomp(t(mat))

percentVar_before <-
  summary(pca_before)$importance[2, 1:2] * 100

df_before <- data.frame(
  PC1 = pca_before$x[,1],
  PC2 = pca_before$x[,2],
  condition = coldata$condition,
  batch = coldata$batch
)

p1 <- ggplot(
  df_before,
  aes(
    PC1,
    PC2,
    color = condition,
    shape = batch
  )
) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = condition),
               linetype = 2) +
  xlab(
    paste0(
      "PC1: ",
      round(percentVar_before[1], 1),
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      round(percentVar_before[2], 1),
      "% variance"
    )
  ) +
  ggtitle("PCA Before SVA") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

print(p1)

# Run SVA
mod  <- model.matrix(~ condition, data=coldata)
mod0 <- model.matrix(~ 1, data=coldata)


svobj <- svaseq(mat, mod, mod0)
svobj$sv
#Correct_batch
mat_sv_corrected <- removeBatchEffect(
  mat,
  covariates = svobj$sv,
  design = mod
)

#PCA_after_SVA
pca_after <- prcomp(t(mat_sv_corrected))

percentVar_after <-
  summary(pca_after)$importance[2, 1:2] * 100

df_after <- data.frame(
  PC1 = pca_after$x[,1],
  PC2 = pca_after$x[,2],
  condition = coldata$condition,
  batch = coldata$batch
)

p2 <- ggplot(
  df_after,
  aes(
    PC1,
    PC2,
    color = condition,
    shape = batch
  )
) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = condition),
               linetype = 2) +
  xlab(
    paste0(
      "PC1: ",
      round(percentVar_after[1], 1),
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      round(percentVar_after[2], 1),
      "% variance"
    )
  ) +
  ggtitle("PCA After SVA") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

print(p2)

#Differential_expression_analysis using sva

coldata$SV1 <- svobj$sv[,1]
coldata$SV2 <- svobj$sv[,2]
coldata$SV3 <- svobj$sv[,3]
coldata$SV4 <- svobj$sv[,4]
coldata$SV5 <- svobj$sv[,5]
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = coldata,
  design = ~ SV1 + SV2 + SV3 + SV4 + SV5 + condition
)
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","B-ALL","Control"))
deg <- as.data.frame(res)

deg_sig <- deg[deg$padj < 0.01 & abs(deg$log2FoldChange) > 1, ]
ensembl_ids <- rownames(deg_sig)
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
deg_sig$symbol <- gene_symbols
#Save_data
write.csv(deg_sig,"SVAseq_DEGs.csv")


