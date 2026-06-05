library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)

#Set thresholds
logFC_cutoff <- 1.5
padj_cutoff  <- 0.01

#Load_Differentially_expressed_DUB_table

deg <- read.xlsx("DEG_DUB.xlsx")


# Classify genes

deg <- deg %>%
  mutate(
    regulation = case_when(
      log2FoldChange >= logFC_cutoff & padj < padj_cutoff ~ "Up",
      log2FoldChange <= -logFC_cutoff & padj < padj_cutoff ~ "Down",
      TRUE ~ "NS"
    )
  )

# Select top genes for labeling

top_genes <- deg %>%
  filter(regulation != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 22)

# Volcano plot
p <- ggplot(deg,
            aes(x = log2FoldChange,
                y = -log10(padj),
                color = regulation)) +
  
  geom_point(
    size = 2,
    alpha = 0.85
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey65",
    linewidth = 0.5
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey65",
    linewidth = 0.5
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = Gene.name,
        color = regulation),
    
    size = 2,
    
    fontface = "plain",
    
    segment.color = NA,
    
    box.padding = 0.2,
    point.padding = 0.15,
    
    max.overlaps = Inf,
    
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Up" = "#2D5A4A",
      "Down" = "#C17C3A",
      "NS" = "grey85"
    )
  ) +
  
  labs(
    title = "Volcano Plot",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~adjusted~italic(P))
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    
    axis.title = element_text(face = "bold"),
    
    legend.position = "top",
    legend.title = element_blank(),
    
    axis.text = element_text(color = "black"),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    )
  )

p

ggsave(
  "Volcano_plot.pdf",
  p,
  width = 4.5,
  height = 5.5,
  dpi = 600
)
