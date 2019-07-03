


library(readxl)
library(sf)
library(raster)
library(dplyr)
library(fasterize)

pts <- read_sf("C:/Users/JG/Desktop/GBIF_DATA/GBIF_Passeriformes_NW_PTES_29N_v1.shp")
vezBasin <- read_sf("D:/MyDocs/R-dev/Thesis/DATA/VECTOR/Bounds/VezBasinCOS12_buff_WGS84_29N.shp")

EFA_MED <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.MED.2001_2016_med.tif")
EFA_DQT <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.DQ1.2001_2016_med.tif")
EFA_SPG <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.SPRG.2001_2016_med.tif")

EFA_MIN <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.Q05.2001_2016_med.tif")
EFA_MAX <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.Q95.2001_2016_med.tif")

EFA_MED[EFA_MED == -3000] <- NA
EFA_DQT[EFA_DQT == 0] <- NA
EFA_SPG[is.na(EFA_MED) | is.na(EFA_DQT)] <- NA
EFA_MIN[EFA_MIN == -3000] <- NA
EFA_MAX[EFA_MAX == -3000] <- NA


#rst <- stack(EFA_MED, EFA_DQT, EFA_SPG, EFA_MIN, EFA_MAX)
#plot(rst)
plot(EFA_MED)
plot(pts[,1],add=TRUE)



for(i in 1:length(rs)){
  
  tmpBuff <- st_buffer(st_as_sf(spObj), dist = rs[i])
  if(i==1){
    all_buffs <- tmpBuff
  }else{
    all_buffs <- rbind(all_buffs,tmpBuff)
  }
}

write_sf(st_cast(all_buffs,"LINESTRING"),"./DATA/GBIF_Buffers_line.shp")
write_sf(all_buffs,"./DATA/GBIF_Buffers_poly.shp")


i <- 0
for(bf_area in c(5,10,15,20)){

  i<-i+1
  distBuff <- rs[bf_area]
  spObj <- SpatialPointsDataFrame(centralCoords, data = data.frame(ID=1:nrow(centralCoords)), 
                                  proj4string = crs("+proj=utm +zone=29 +datum=WGS84 +units=m +no_defs"))
  buffs <- st_buffer(st_as_sf(spObj), dist = distBuff)
  
  rstBuffs <- fasterize(buffs,raster = EFA_MED, field = "ID")
  
  rst <- stack(rstBuffs, EFA_MED, EFA_DQT, EFA_SPG, EFA_MIN, EFA_MAX)
  
  rstDF <- values(rst) %>% na.omit
  colnames(rstDF) <- c("ID", "EFA_MED", "EFA_DQT", "EFA_SPG", "EFA_MIN", "EFA_MAX")
  
  rstDF_agg <- rstDF %>% 
    as.data.frame %>% 
    group_by(ID) %>% 
    summarise_all(.funs = list("STD"=sd), na.rm=TRUE) %>% 
    arrange(ID) %>% 
    cbind(spRichSAR %>% filter(a == bf_area) %>% select(s_avg))
    # cbind(spRichSAR %>% filter(a == 5) %>% select(s_avg) %>%  rename("s_avg_5k"  = "s_avg")) %>% 
    # cbind(spRichSAR %>% filter(a == 10) %>% select(s_avg) %>% rename("s_avg_10k" = "s_avg")) %>% 
    # cbind(spRichSAR %>% filter(a == 15) %>% select(s_avg) %>% rename("s_avg_15k" = "s_avg")) %>% 
    # cbind(spRichSAR %>% filter(a == 20) %>% select(s_avg) %>% rename("s_avg_20k" = "s_avg"))
  
  #cor(rstDF_agg %>% select(s_avg), rstDF_agg %>% select(2:7)) %>% round(2)
  
  tmpCor <- cor(rstDF_agg %>% 
        select(s_avg), rstDF_agg %>% 
        select(2:6), method="spearman") %>% 
    round(2)
  
  tmpCor <- data.frame(a=bf_area, tmpCor)
  print(tmpCor)
  
  if(i==1){
    corMatSpear <- tmpCor
  }else{
    corMatSpear <- rbind(corMatSpear, tmpCor)
  }
  
}



