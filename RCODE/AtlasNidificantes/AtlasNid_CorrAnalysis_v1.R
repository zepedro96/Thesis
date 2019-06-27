

library(readxl)
library(sf)
library(raster)
library(dplyr)
library(fasterize)

atlasDF <- read_excel("D:/MyDocs/Dropbox/MasterThesisJoséSilva/DATA/AtlasAvesNidificantes/Avifauna_AlllRecordsLabeled_UTM10k.xls")

grid10k <- read_sf("D:/MyDocs/Dropbox/MasterThesisJoséSilva/DATA/AtlasAvesNidificantes/UTMgrid10k_projSIMBioN/UTMgrid10k_UTM29N_pSIMB_v1.shp")


countDistinct <- function(x) length(unique(x))


atlasDF_spR <- atlasDF %>% 
                group_by(utm_10k) %>% 
                summarize(spRich = countDistinct(spabrev))


grid10k <- grid10k %>% full_join(atlasDF_spR, by=c("mgrs10K"="utm_10k"))

write_sf(grid10k,"./OUT/Grid10k_SpRich.shp")

plot(grid10k)

EFA_MED <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.MED.2001_2016_med.tif")
EFA_DQT <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.DQ1.2001_2016_med.tif")
EFA_SPG <- raster("D:/MyDocs/R-dev/PA_EcosystemStability/DATA/RASTER/MODIS/PT/IND/MultiAnnual/MOD13Q1.EVI.SPRG.2001_2016_med.tif")

gridRaster <- fasterize(grid10k, EFA_MED, field = "OBJECTID")
crs(gridRaster) <- crs(EFA_MED)

rstDF <- stack(EFA_MED, EFA_DQT, EFA_SPG, gridRaster) %>% values %>% na.omit
colnames(rstDF) <- c("MED","DQT","SPG","OBJECTID")

rstDF_Agg10k <- rstDF %>% 
  as.data.frame %>% 
  group_by(OBJECTID) %>% 
  summarize_all(.funs=list("STD"=sd, "AVG"=mean))

rstDF_All <- rstDF_Agg10k %>% 
  inner_join(grid10k, by="OBJECTID") %>% 
  st_as_sf %>% 
  # st_drop_geometry() %>% 
  select(-area_m2, -Shape_Leng, -Shape_Area)


plot(rstDF_All)

plot(rstDF_All %>% select(2:7,9) %>% st_drop_geometry())

cor(rstDF_All %>% select(2:7,9) %>% st_drop_geometry() %>% na.omit, method="spearman") %>% round(2)

