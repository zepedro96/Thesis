

library(raster)
library(fasterize)
library(dplyr)
library(trend)
library(readxl)
library(sf)
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



senSlope <- function(x,...) sens.slope(x,...)$estimates
senSlopePval <- function(x,...) sens.slope(x,...)$p.value

grid1km   <- read_sf("D:/Dropbox/MasterThesisJoséSilva/DATA/Grids/SampGridVez_1km_Vez_All_v1.shp")

## Continuous EFA indicators
## EFA_max / EFA_semax

EFA_rst <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAavg.tif")

EFA_name <- "av" ## USE TWO LETTERS ONLY (shapefiles will not work with more characters... :-0

names(EFA_rst) <- paste(EFA_name,substr(as.character(1984:2019),3,4),sep="_")

birdGridRst1km <- fasterize::fasterize(grid1km, EFA_rst[[1]], field = "ID_PSU")

rstStack <- stack(birdGridRst1km, EFA_rst)

rstDF <- values(rstStack) %>% na.omit
colnames(rstDF)[1] <- "ID_PSU"
  
rstDF_agg <- rstDF %>% 
  as.data.frame %>% 
  group_by(ID_PSU) %>% 
  summarize_all(.funs=list("sd"=sd,"av"=mean)) 



## Estimate trends over different periods -------------------------------------------- ##

# Last ten years trend
trend_10yr <- rstDF_agg %>% 
  select(27:36) %>% 
  data.frame(sen_10yr = apply(.,1,senSlope),pv_10yr = apply(.,1,senSlopePval)) %>% 
  select(sen_10yr, pv_10yr)

trend_20yr <- rstDF_agg %>% 
  select(17:36) %>% 
  data.frame(sen_20yr = apply(.,1,senSlope),pv_20yr = apply(.,1,senSlopePval)) %>% 
  select(sen_20yr, pv_20yr)

trend_all <- rstDF_agg %>% 
  select(3:36) %>% 
  data.frame(sen_all = apply(.,1,senSlope),pv_all = apply(.,1,senSlopePval)) %>% 
  select(sen_all, pv_all)


trendsDF <- data.frame(ID_PSU=rstDF_agg$ID_PSU, trend_10yr, trend_20yr, trend_all)

grid1km <- grid1km %>% 
  left_join(rstDF_agg, by = "ID_PSU")%>% 
  left_join(trendsDF, by = "ID_PSU")


write_sf(grid1km,"./OUTtoSHARE/EFA_trends/EFA_avg_grid1km_1984_2019_v1.shp")

write_sf(grid1km %>% st_centroid,"./OUTtoSHARE/EFA_trends/EFA_avg_grid1kmCentroids_1984_2019_v1.shp")


## Calculate Shannon

