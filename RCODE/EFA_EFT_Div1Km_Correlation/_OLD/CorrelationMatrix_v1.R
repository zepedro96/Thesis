
vezDF <- readRDS("./DATA/vezDF_GHC_merge.rds")

selCols <- c(51:56,4,19,20,21,41:49)
#selColNames<- c("spRich")

cmp <- round(cor(vezDF[,selCols]), 2)

cms <- round(cor(vezDF[,selCols], method="spearman"), 2)

cmp[1:9,]

cms[1:9,]
