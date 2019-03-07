

#removing burnt plots from data frame

burnt_plots <- read.csv("D:/Thesis/DATA/fire_birdsgrid.csv")

tbl_data <- read.csv("D:/Thesis/DATAtoShare/winter15_FNL.csv")

dlt_burnt <- tbl_data[ !(tbl_data$ID_SSU %in% burnt_plots$ID_SSU), ]

write.csv(er_burnt, "./OUT/dlt_burnt_winter15.csv")
