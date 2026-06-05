library(tidyverse)
library(broom)
library(ggplot2)
library(survminer)
library(survival)

#Upload_Clinical_file_having_all_clinical_parameters

clin <- read.csv("clinical_data.csv")
clin <- clin[,-1]

#Set_lowOTUB1_expression_group_as_reference

clin$Group<- factor(
  clin$Group,
  levels = c("Low", "High")
)
#Perform_multivariate_Cox_regression
cox <- coxph(Surv(Time, Status) ~ Group + age + Gender + Molecular_Subtype + MRD...at.day.29 + WBC,
              data = clin)
ggforest(cox, data = clin)

#Convert_into_tidy_dataframe
df <- tidy(cox, exponentiate = TRUE, conf.int = TRUE)
df$term <- c("OTUB1 expression", "Age", "Gender", "BCR-ABL1", "ETV6-RUNX1", "Hyperdiploidy", "iamp21", "MLL-rearranged", "TCF3 fusion", "MRD @ day 29", "WBC")
df$term <- factor(df$term, levels = rev(df$term))

df$var <- df$term

custom_colors <- c("OTUB1 expression" ="#D95F02", "Age" ="#B8860B", "Gender" ="#009E73", "BCR-ABL1" ="#B22222", "ETV6-RUNX1" ="#B22222", "Hyperdiploidy" ="#B22222", "iamp21" ="#B22222", "MLL-rearranged" ="#B22222", "TCF3 fusion" ="#B22222" , "MRD @ day 29" ="#6A3D90", "WBC" ="#C21857")

#Create_label_text(HR + CI + p-value)

df$label <- sprintf(
  "%.2f (%.2f–%.2f)\np %s",
  df$estimate,
  df$conf.low,
  df$conf.high,
  format.pval(df$p.value, digits = 2, eps = 0.001)
)
#Make_final_forestplot
ggplot(df, aes(x = estimate, y = term, fill = var, color = var)) +
  
  # Confidence_intervals
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, size = 0.8) +
  
  #Points
  geom_point(size = 3.5, shape = 21, color = "black") +
  
  #Reference_line_at_HR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "#654321", size = 0.8) +
  
  # Text_labels_on_the_right
  geom_text(aes(label = label),
            hjust = -0.1, size = 4, family = "sans") +
  
  # Log_scale
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50)
  ) +
  coord_cartesian(xlim = c(0.25, 4), clip = "off")+
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
 
   # Titles_and_labels
  
  labs(
    x = "Hazard Ratio (log scale)",
    y = NULL,
    title = "Multivariable Cox Regression",
    subtitle = "Hazard ratios with 95% confidence intervals"
  ) +

  coord_cartesian(clip = "off") +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    plot.margin = margin(5.5, 80, 5.5, 5.5),  # extra space for labels
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
