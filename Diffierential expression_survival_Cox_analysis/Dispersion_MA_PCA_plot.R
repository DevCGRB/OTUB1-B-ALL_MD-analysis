library(ggplot2)
library(ggrastr)
library(patchwork)

 #Create_common_theme_for_multipanel_figure

common_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_blank(),
    
    axis.title = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.text = element_text(size = 10),
    
    legend.title = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )

#Dispersion_Plot

disp_df <- data.frame(
  mean      = mcols(dds)$baseMean,
  gene_est  = mcols(dds)$dispGeneEst,
  final_est = mcols(dds)$dispersion,
  fitted    = mcols(dds)$dispFit
)

disp_df <- disp_df[
  complete.cases(disp_df) &
    disp_df$mean > 0 &
    disp_df$gene_est > 0 &
    disp_df$final_est > 0 &
    disp_df$fitted > 0,
]

disp_df <- disp_df[order(disp_df$mean), ]

Disp_P <- ggplot(disp_df, aes(mean)) +
  
  rasterise(
    geom_point(
      aes(
        y = gene_est,
        color = "Gene estimate"
      ),
      alpha = 0.20,
      size = 0.5
    ),
    dpi = 600
  ) +
  
  rasterise(
    geom_point(
      aes(
        y = final_est,
        color = "Final estimate"
      ),
      alpha = 0.35,
      size = 0.6
    ),
    dpi = 600
  ) +
  
  geom_line(
    aes(
      y = fitted,
      color = "Fitted"
    ),
    linewidth = 0.8
  ) +
  
  scale_x_log10() +
  scale_y_log10() +
  
  scale_color_manual(
    values = c(
      "Gene estimate" = "#4DAF4A",
      "Final estimate" = "#D95F02",
      "Fitted" = "#7570B3"
    )
  ) +
  
  labs(
    x = "Mean normalized counts",
    y = "Dispersion"
  ) +
  
  common_theme +
  
  theme(
    legend.position = "none"
  )

#MA_Plot

res <- results(
  dds,
  contrast = c(
    "condition",
    "B-ALL",
    "controls"
  )
)

res <- as.data.frame(res)

res <- res[
  !is.na(res$padj) &
    !is.na(res$baseMean) &
    is.finite(res$log2FoldChange) &
    res$baseMean > 0,
]

res_plot <- subset(
  res,
  baseMean > 50
)

res_plot$significant <- factor(
  ifelse(
    res_plot$padj < 0.01 &
      abs(res_plot$log2FoldChange) > 1.5,
    "Significant",
    "Not significant"
  )
)


MA_P <- ggplot(
  res_plot,
  aes(
    baseMean,
    log2FoldChange
  )
) +
  
  rasterise(
    geom_point(
      aes(color = significant),
      alpha = 0.5,
      size = 0.6
    ),
    dpi = 600
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  scale_x_log10() +
  
  scale_color_manual(
    values = c(
      "Not significant" = "grey75",
      "Significant" = "#D73027"
    )
  ) +
  
  coord_cartesian(
    ylim = c(-5, 5)
  ) +
  
  labs(
    x = "Mean normalized counts",
    y = expression(log[2]~Fold~Change)
  ) +
  
  common_theme +
  
  theme(
    legend.position = "none"
  )

#PCAplot_afternormalization

vsd <- vst(dds, blind = TRUE)

mat <- assay(vsd)

pca <- prcomp(t(mat))

percentVar <- summary(pca)$importance[
  2,
  1:2
] * 100

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  condition = colData(dds)$condition
)

PCA_P <- ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    color = condition
  )
) +
  
  geom_point(
    size = 3.5,
    alpha = 0.9
  ) +
  
  stat_ellipse(
    aes(group = condition),
    linewidth = 0.7,
    linetype = 2
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(percentVar[1],1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(percentVar[2],1),
      "%)"
    )
  ) +
  
  common_theme +
  
  theme(
    legend.position = "bottom"
  )
FinalFig <- (
  Disp_P + MA_P
) /
  PCA_P +
  plot_layout(
    heights = c(1, 1.4)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 14
      )
    )
  )

FinalFig

ggsave(
  "FinalFig.pdf",
  FinalFig,
  width = 8,
  height = 7,
  units = "in",
  device = cairo_pdf
)
