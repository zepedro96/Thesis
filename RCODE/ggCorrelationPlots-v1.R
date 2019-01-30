

library(ggcorrplot)
library(ggplot2)


vezDF <- readRDS("./DATAtoShare/vezDF_GHC_merge.rds")


selCols <- c(51:56,4,19,20,21,41:49)
#selColNames<- c("spRich")

cmp <- round(cor(vezDF[,selCols]), 2)

cms <- round(cor(vezDF[,selCols], method="spearman"), 2)

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(vezDF[,selCols])

g <- ggcorrplot(cms, #hc.order = TRUE,
           type = "lower", #p.mat = p.mat, 
           lab = TRUE,
           digits = 1)

ggsave(filename="./OUT/corrMatrix_Spring16.png", plot = g, width = 8, height=8)
