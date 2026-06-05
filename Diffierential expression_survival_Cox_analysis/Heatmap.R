library(openxlsx)
library(circlize)
library(ComplexHeatmap)

#Input_differentialDUB_data_having__expression_ofeach_samples_in_columns

df <- read.xlsx("DUB expression.xlsx")

rownames(df) <- make.names(df$Gene.name, unique = T)

df <- df[,-1]

#Convert each DUB's expression values into row-wise z-scores

mat1.z <- t(apply(df, 1, scale))


#Set_the_color_ramp

col_fun <- colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027"))

#Draw_Heatmap

h1 <- Heatmap(mat1.z, name = "mat", col = col_fun)
draw(h1)
