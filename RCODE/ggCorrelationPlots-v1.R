

library(ggcorrplot)
library(ggplot2)


vezDF <- read.csv("./DATA/dlt_burnt_kmeans_ndvi16.csv")

vezDF <-vezDF[,-c(1)]

selCols <- c(42:48,4,19,20,21,40,41)
#selColNames<- c("spRich")

cmp <- round(cor(vezDF[,selCols]), 2)

cms <- round(cor(vezDF[,selCols], method="spearman"), 2)

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(vezDF[,selCols])

g <- ggcorrplot(cms, #hc.order = TRUE,
           type = "lower", #p.mat = p.mat, 
           lab = TRUE,
           digits = 1)

ggsave(filename="./OUT/corrMatrix_k_ndvi16.png", plot = g, width = 8, height=8)
