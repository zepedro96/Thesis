

library(readxl)
library(dplyr)
library(cluster)

source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r")

x <- read_excel("C:/Users/JG/Downloads/R_sheet_birds_v2.xlsx")

xc <- cor(x[,c(26, 28, 30)]) %>% round(2)

km <- kmeans(x[,c(26,28)], 4)

boxplot(x$Size_avg)
boxplot(x$Weight_avg)

pamClust <- pam(x,k=4,metric = "euclidean",diss=FALSE)

km$cluster

boxplot.with.outlier.label(x$Size_avg, x$Species_name, main="Average size (cm)", ylab="Size (cm)")
boxplot.with.outlier.label(x$Weight_avg, x$Species_name, main="Average weight (g)", ylab="Weight (g)")


