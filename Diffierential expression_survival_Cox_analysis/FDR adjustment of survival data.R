library(openxlsx)

#input_logrankp_data_from_survival_genie

df <- read.xlsx("survival_logrankp.xlsx")

df <- as.data.frame(df)

#Adjust_FDR

df$peripheral_FDR <- p.adjust(df$LogRankPB, method = "BH")
df$bonemarrow_FDR <- p.adjust(df$LogRankBM, method = "BH")

df$pb_significant <- df$peripheral_FDR < 0.05
df$bm_significant <- df$bonemarrow_FDR < 0.05

df
