

#removing burnt plots from data frame

burnt_plots <- read.csv("D:/Thesis/DATA/fire_birdsgrid.csv")

tbl_data <- read.csv("D:/Thesis/DATAtoShare/TABLES/autumn16_FNL.csv")

dlt_burnt <- tbl_data[ !(tbl_data$ID_SSU %in% burnt_plots$ID_SSU), ]

write.csv(dlt_burnt, "./DATAtoShare/TABLES/dlt_burnt_autumn16.csv")

