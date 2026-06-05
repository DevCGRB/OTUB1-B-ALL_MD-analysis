library(DESeq2)
library(readr)
library(dplyr)
library(survival)
library(survminer)
library(broom)

#Extract_OTUB1/ENSG00000167770_expression_from_raw_TARGET-ALL-P2_counts

counts <- read.csv("TARGETP2count.csv", row.names = 1)
colnames(counts) <- gsub("\\.", "-", colnames(counts))

colData <- data.frame(
  row.names = colnames(counts),
  condition = rep("B_all", ncol(counts))  # placeholder
)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ 1)

vsd <- vst(dds, blind = TRUE)

norm_mat <- assay(vsd)
gene_id <- "ENSG00000167770"

gene_expr <- norm_mat[gene_id, ]

expr_df <- data.frame(
  ID = names(gene_expr),
  gene_expr = as.numeric(gene_expr)
)

#Merge_OTUB1expression_in_clinicaldata_file

clinical <- read.csv("clinical information.csv")
clinical$time <- ifelse(clinical$vital_status == "Dead",
                        clinical$days_to_death,
                        clinical$days_to_last_follow_up)
clinical$status <- ifelse(clinical$vital_status == "Dead", 1, 0)

clinical <- clinical[!is.na(clinical$time), ]
clinical$age <- clinical$age_at_diagnosis / 365

get_patient_id <- function(x) {
  sapply(strsplit(x, "-"), function(y) paste(y[1:3], collapse = "-"))
}
merged <- merge(clinical, expr_df, by = "ID")
merged$group <- ifelse(merged$gene_expr > median(merged$gene_expr),
                       "High", "Low")
merged$patient_id <- get_patient_id(merged$ID)

#Perform_Univariate_Cox_analysis

cox1 <- coxph(Surv(time, status) ~ gene_expr, data = merged)
summary(cox1)

ggforest(cox1, data = merged)

#Extract_results

cox.res <- tidy(
  cox1,
  exponentiate = TRUE,
  conf.int = TRUE
)

#Create_labels

cox.res$HR_label <- paste0(
  round(cox.res$estimate, 2),
  " (",
  round(cox.res$conf.low, 2),
  "-",
  round(cox.res$conf.high, 2),
  ")"
)

#Create_forest_plot

p <- ggplot(cox.res,
            aes(x = term,
                y = estimate,
                ymin = conf.low,
                ymax = conf.high)) +
  
  geom_pointrange(size = 1.2) +
  
  geom_hline(yintercept = 1,
             linetype = "dashed",
             linewidth = 0.7) +
  
  coord_flip() +
  
  geom_text(aes(label = HR_label,
                y = conf.high + 0.3),
            size = 5) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "",
    y = "Hazard Ratio (95% CI)",
    title = "Univariate Cox Regression Analysis"
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_blank()
  )

print(p)



