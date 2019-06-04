
library(raster)
library(sf)
library(fasterize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RStoolbox)
library(psych)
library(microbiome)


divEvFun <- function(x){
  as.data.frame(
    c(evenness(as.numeric(table(x)), index = "all", zeroes = FALSE, detection = 0)[1,,drop=TRUE],
      diversity(as.numeric(table(x)), index = c("shannon","inverse_simpson","gini_simpson"), zeroes = FALSE)[1,,drop=TRUE],
      eft_count = length(unique(x)))
  )
}

shannon <- function (x, na.rm=TRUE){
  x <- na.omit(x)
  if(length(x)==0){
    return(0)
  }else{
    tb <- table(x)
    tbrel <- tb / sum(tb)
    z <- tbrel * log(tbrel)
    out <- -1 * sum(z)
    return(out)
  }
}

countDistinct <- function(x){
  x <- na.omit(x)
  if(length(x)==0){
    return(0)
  }else{
    return(length(unique(x))) 
  }
} 

makePivotTableRelFreq <- function(x,ID_field,Data_field){
  
  pivotTable <- x %>% 
    group_by(.data[[ID_field]], .data[[Data_field]]) %>% 
    summarize(nc = n()) %>% 
    spread(key = Data_field, value = "nc", fill = 0)
  
  out <- data.frame(
    pivotTable[,ID_field, drop=FALSE],
    pivotTable[,-1] / apply(pivotTable[,-1],1,sum))
  
  return(out)
}

getVarNames <- function(x) gsub(".tif","",unlist(lapply(strsplit(fpaths,"_"),function(x) x[length(x)])))


birdsgrid <- read_sf("D:/MyDocs/Dropbox/MasterThesisJoséSilva/DATA/FieldData/Birds/Grids/Grid200mEffectSampSSU_WGS84UTM29N_VEZ_v1.shp")
grid1km   <- read_sf("D:/MyDocs/GeoData/ByProject/IND_CHANGE/FieldWorkData/CommonBirds_FieldData_MainTable/Grids/Grid1kmEffectSampPSU_WGS84UTM29N_VEZ_v1.shp")

fpaths <- list.files("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun", pattern=".tif$", full.names = TRUE)[-c(13:14)]

vnames <- getVarNames(fpaths)

for(i in 1:length(fpaths)){
  
  print(i)
  tmp <- stack(fpaths[[i]])[[31]] # yr=31 -> 2014
  if(i==1)
    rst <- stack(tmp)
  else
    rst <- stack(rst, tmp)
}

names(rst) <- vnames

# EFT map / discrete combination of EFA
map2 <- raster("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFTcomnbs01.tif")



plot(map2)

## ------------------------------------------------------------------------------------- ##

birdGridRst <- fasterize::fasterize(birdsgrid, rst[[1]], field = "ID_SSU")

birdGridRst1km <- fasterize::fasterize(grid1km, rst[[1]], field = "ID_PSU")

## ------------------------------------------------------------------------------------- ##

rstDF <- na.omit(values(stack(birdGridRst, rst))) %>% 
  as.data.frame %>% 
  rename("ID_SSU" = "layer") %>% 
  group_by(ID_SSU) %>% 
  summarize_all(.funs=list(avg = mean, std = sd))

rstDF_km20c <- na.omit(values(stack(birdGridRst, map2))) %>% 
  as.data.frame %>% 
  rename("ID_SSU" = "layer") %>% 
  group_by(ID_SSU) %>% 
  #summarize_all(.funs=list(shn20c = shannon, nclass20c = countDistinct))
  do(divEvFun(.[[2]])) %>% 
  ungroup %>% 
  select(-1) %>% 
  cor %>% round(2)


DFkm20c <- na.omit(values(stack(birdGridRst, map2))) %>% 
  as.data.frame %>% 
  `colnames<-`(c("ID_SSU","clust")) %>% 
  makePivotTableRelFreq("ID_SSU","clust")

rstDF <- rstDF %>% 
  left_join(rstDF_km20c, by="ID_SSU") %>% 
  left_join(DFkm20c, by="ID_SSU")

## 1km data Primary Sampling Units (PSU) --------

rstDF1km <- na.omit(values(stack(birdGridRst1km, rst))) %>% 
  as.data.frame %>% 
  rename("ID_PSU" = "layer") %>% 
  group_by(ID_PSU) %>% 
  summarize_all(.funs=list(avg = mean, std = sd))

rstDF_km20c_1km <- na.omit(values(stack(birdGridRst1km, map2))) %>% 
  as.data.frame %>% 
  rename("ID_PSU" = "layer") %>% 
  group_by(ID_PSU) %>% 
  #summarize_all(.funs=list(shn20c = shannon, nclass20c = countDistinct))
  do(divEvFun(.[[2]]))

DFkm20c_1km <- na.omit(values(stack(birdGridRst1km, map2))) %>% 
  as.data.frame %>% 
  `colnames<-`(c("ID_PSU","clust")) %>% 
  makePivotTableRelFreq("ID_PSU","clust")

rstDF1km <- rstDF1km %>% 
  left_join(rstDF_km20c_1km, by="ID_PSU") %>% 
  left_join(DFkm20c_1km, by="ID_PSU")


## ------------------------------------------------------------------------------------- ##

vezDF <- readRDS("./DATAtoShare/RDATA/vezDF_GHC_merge.rds")

vezDF1km <- read.csv("D:/MyDocs/GeoData/ByProject/IND_CHANGE/FieldWorkData/CommonBirds_FieldData_MainTable/birdDiv_multiSARavg_preds_Vez_v2.csv")

vezDF_vars <- vezDF %>% 
  left_join(rstDF, by=c("ID_SSU.x"="ID_SSU"))

vezDF_vars1km <- vezDF1km %>% 
  left_join(rstDF1km, by=c("ID_PSU"))

cmPears <- cor(vezDF_vars[,-c(1:3)]) %>% round(2)
cmSpear <- cor(vezDF_vars[,-c(1:3)], method="spearman") %>% round(2)

cmPears_1km <- cor(vezDF_vars1km[,-1]) %>% round(2)
cmSpear_1km <- cor(vezDF_vars1km[,-1], method="spearman") %>% round(2)




write.csv(cmPears, "./OUT/cmPears_LT-vars-jul-jun-v4.csv")
write.csv(cmSpear, "./OUT/cmSpear_LT-vars-jul-jun-v4.csv")
write.csv(cmPears_1km, "./OUT/cmPears_1km_LT-vars-jul-jun-v4.csv")
write.csv(cmSpear_1km, "./OUT/cmSpear_1km_LT-vars-jul-jun-v4.csv")


## ------------------------------------------------------------------------- ##



g1 <- ggplot(vezDF_vars, aes(x=EFAamp_std, y=spRich)) + 
  geom_point() + 
  geom_smooth(method="lm")
  #geom_quantile(quantiles=c(0.1,0.5,0.9))

plot(g1)

g2 <- ggplot(vezDF_vars, aes(x=ESPI_std, y=spRich)) + 
  geom_point() + 
  #geom_smooth(method="lm") + 
  geom_quantile(quantiles=seq(0.1,0.9,0.1))

plot(g2)


## ------------------------------------------------------------------------- ##





