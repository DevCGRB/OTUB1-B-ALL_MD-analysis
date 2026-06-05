library(DESeq2)
library(ggrastr)
library(dplyr)
library(patchwork)
library(AnnotationDbi)
library(org.Hs.eg.db)

#Load_combined_Rawcounts_file
count_data <- read.csv("TARGETP2vsControl_merged.csv", row.names = 1)
#remove_low_count_genes
count_data <- count_data[which(rowSums(count_data) > 100),]
#Segregate_Disease_and_Control_samples
condition <- factor(c(rep("B-ALL", 198), rep("controls", 24)))

coldata <- data.frame(row.names = colnames(count_data), condition)

#Create_DESeq2_object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design=~condition)


#Deseq2 analysis
dds <- DESeq(dds)

#Fetch Deseq2 results

res <- results(dds, contrast = c("condition", "B-ALL", "controls"))

res <- res[!is.na(res$padj) & is.finite(res$log2FoldChange), ]

#Save differential expression results
res <- as.data.frame(res)
#Filter only significant genes
resdatafilt <- res[(res$padj < 0.01) & (res$baseMean > 50) & (abs(res$log2FoldChange) > 1.5),]

resdatafilt <- resdatafilt[order(resdatafilt$log2FoldChange, decreasing = TRUE),]
#Convert ensembl IDs into gene names
ensembl_ids <- rownames(resdatafilt)
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
resdatafilt$geneid <- gene_symbols
#save_result
write.csv(resdatafilt,"DEGresult.csv")

