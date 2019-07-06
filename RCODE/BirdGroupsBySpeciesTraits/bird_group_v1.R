library(readxl)
library(lattice)
#Birds grouping

table_bd <-read_xlsx("./DOCS/R_sheet_birds.xlsx")

qts <- quantile(table_bd$Size_avg, na.rm=TRUE)

ct <- cut(table_bd$Size_avg, qts, labels = 1:4, include.lowest = TRUE, include.highest = TRUE)


          