

library(readxl)
library(sf)
library(raster)
library(dplyr)
library(fasterize)

atlasDF <- read_excel("D:/MyDocs/Dropbox/MasterThesisJoséSilva/DATA/AtlasAvesNidificantes/Avifauna_AlllRecordsLabeled_UTM10k.xls")

grid10k <- read_sf("D:/MyDocs/Dropbox/MasterThesisJoséSilva/DATA/AtlasAvesNidificantes/UTMgrid10k_projSIMBioN/UTMgrid10k_UTM29N_pSIMB_v1.shp")

vezBasin <- read_sf("D:/MyDocs/R-dev/Thesis/DATA/VECTOR/Bounds/VezBasinCOS12_buff_WGS84_29N.shp")

st_crs(grid10k) <- st_crs(vezBasin)

grid10kVez <- grid10k %>% filter(st_intersects(grid10k, st_buffer(vezBasin, dist = 5000), 
                                               sparse = FALSE)[,1])


countDistinct <- function(x) length(unique(x))


atlasDF_spR <- atlasDF %>% 
                group_by(utm_10k) %>% 
                summarize(spRich = countDistinct(spabrev))


grid10kVez <- grid10kVez %>% full_join(atlasDF_spR, by=c("mgrs10K"="utm_10k"))

#write_sf(grid10kVez,"./OUT/Grid10k_SpRich.shp")

plot(grid10kVez)

EFA_MED <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.MED.2001_2016_med.tif")
EFA_DQT <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.DQ1.2001_2016_med.tif")
EFA_SPG <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.SPRG.2001_2016_med.tif")

gridRaster <- fasterize(grid10kVez, EFA_MED, field = "OBJECTID")
crs(gridRaster) <- crs(EFA_MED)

rstDF <- stack(EFA_MED, EFA_DQT, EFA_SPG, gridRaster) %>% values %>% na.omit
colnames(rstDF) <- c("MED","DQT","SPG","OBJECTID")

rstDF_Agg10k <- rstDF %>% 
  as.data.frame %>% 
  group_by(OBJECTID) %>% 
  summarize_all(.funs=list("STD"=sd, "AVG"=mean))

rstDF_All <- rstDF_Agg10k %>% 
  inner_join(grid10kVez, by="OBJECTID") %>% 
  st_as_sf %>% 
  # st_drop_geometry() %>% 
  select(-area_m2, -Shape_Leng, -Shape_Area) #%>% 
  #filter(!(OBJECTID %in% c(941, 963)))

#write_sf(rstDF_All,"./OUT/Grid10k_SpRich_VezOnly_v2.shp")



plot(rstDF_All %>% select(-OBJECTID, -mgrs10K))

plot(rstDF_All %>% select(2:7,9) %>% st_drop_geometry())

cor(rstDF_All %>% select(2:7,9) %>% 
      st_drop_geometry() %>% na.omit, method="spearman") %>% round(2)

plot(rstDF_All$DQT_STD, rstDF_All$spRich)


