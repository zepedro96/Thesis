
library(raster)
library(fasterize)
library(sf)
library(dplyr)

vez <- read_sf("./DATA/VECTOR/Bounds/VezBasinCOS12_WGS84_29N.shp")

r <- raster("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/EFT/LT578comp_EVI_jul-jun_EFT_2015_med_amp_semax_025qts_v1.tif")

vezRst <- fasterize(vez, r)

rMsk <- mask(r, vezRst)

v <- na.omit(values(rMsk))
tb <- table(v)

tb <- sort(tb, decreasing=TRUE)

length(unique(v))

EFT_DF <- data.frame(EFT_class = names(tb),
           EFT_Prod = substr(names(tb),1,1),
           EFT_Seas = substr(names(tb),2,2),
           EFT_Phen = substr(names(tb),3,3),
           PercCover = (as.numeric(tb)/sum(tb))*100)

EFT_DF1 <- EFT_DF %>% group_by(EFT_Phen) %>% summarize(EFT_Prod_sum = sum(PercCover))


write.csv(EFT_DF, "./OUTtoShare/EFT_PercFrequencyTable-v1.csv", row.names = FALSE)
write.csv(EFT_DF1, "./OUTtoShare/EFT_Phenology_PercFrequencyTable-v1.csv", row.names = FALSE)

