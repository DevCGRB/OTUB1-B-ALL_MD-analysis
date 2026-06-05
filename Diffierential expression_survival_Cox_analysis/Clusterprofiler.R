library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)


#Input_listof_Differentially_expressed_DUBs

df <- read.xlsx("DEG_DUBs.xlsx")


#Recheck_significance

sig_genes <- df$Gene.name[
  df$adj.p.value < 0.01 &
    abs(df$log2FoldChange) > 1.5
]

#Convert_gene_symbols_ito_ENTREZid

gene_map <- bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

entrez_ids <- unique(gene_map$ENTREZID)

cat("Mapped genes:", length(entrez_ids), "\n")


#Reactome_enrichment

reactome_res <- enrichPathway(
  gene = entrez_ids,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.20,
  readable = TRUE
)

head(as.data.frame(reactome_res))

#Reactome_Dot_Plot

p <- dotplot(
  reactome_res,
  showCategory = 10,
  x = "GeneRatio",
  color = "p.adjust"
) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    )
  )

print(p)
# Save figure
ggsave(
  "Reactome_Dotplot.pdf",
  p,
  width = 10,
  height = 8
)

