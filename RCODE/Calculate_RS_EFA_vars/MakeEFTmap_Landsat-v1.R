
library(raster)
library(sf)
library(fasterize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RStoolbox)
library(microbiome)

getVarNames <- function(x) gsub(".tif","",unlist(lapply(strsplit(fpaths,"_"),function(x) x[length(x)])))

fpaths <- list.files("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun", pattern=".tif$", full.names = TRUE)[]


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

## ------------------------------------------------------------------------------------------------------- ##

avg <- values(rst[["EFAmed"]])
std <- values(rst[["EFAamp"]])
smx <- values(rst[["EFAsemax"]])

cor(cbind(avg,std,smx), method="spearman") %>% round(2)

qts <- seq(0, 1, by = 0.25)

avgQt <- quantile(avg, qts, na.rm = TRUE)
stdQt <- quantile(std, qts, na.rm = TRUE)

recDF <- data.frame(avgCl = as.character(cut(avg, avgQt, labels = 1:4)),
                    stdCl = as.character(cut(std, stdQt, labels = 1:4)),
                    mmax  = as.character(smx))

funConcat <- function(x, na.rm=TRUE){
  if(any(is.na(x))){
    return(NA)
  }else{
    return(as.integer(paste(x,collapse="",sep="")))
  }
}

EFT_combns <- apply(recDF, MARGIN = 1, FUN = funConcat)

length(unique(EFT_combns))

EFT_map <- raster(rst[[1]])
values(EFT_map) <- EFT_combns

plot(EFT_map)

writeRaster(EFT_map,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFTcomnbs.tif", 
            overwrite = TRUE)



