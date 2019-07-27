
library(raster)
library(sf)
library(fasterize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RStoolbox)
library(microbiome)

funConcat <- function(x, na.rm=TRUE){
  if(any(is.na(x))){
    return(NA)
  }else{
    return(as.integer(paste(x,collapse="",sep="")))
  }
}

getVarNames <- function(x) gsub(".tif","",unlist(lapply(strsplit(fpaths,"_"),function(x) x[length(x)])))



fpaths <- list.files("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun", pattern=".tif$", full.names = TRUE)

vnames <- getVarNames(fpaths)


# Output folder
outFolder <- "./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/EFT"
# Use quartiles
qts <- seq(0, 1, by = 0.25)
# Output years
yrs <- 1984:2019


## ---------------------------------------------------------------- ##

for(in_yr in 1:length(yrs)){
  
  yr <- yrs[in_yr]
  
  # Load each EFA in the input folder by year 
  #
  pb <- txtProgressBar(min=1,max=length(fpaths),style=3)
  
  for(i in 1:length(fpaths)){
    
    # Load one EFA variable and select the i-th year
    #
    tmp <- stack(fpaths[i])[[in_yr]] # yr=31 -> 2014
    
    if(i==1)
      rst <- stack(tmp)
    else
      rst <- stack(rst, tmp)
    
    setTxtProgressBar(pb, i)
  }
  
  names(rst) <- vnames
  
  ## ---------------------------------------------------------------- ##
  
  avg <- values(rst[["EFAmed"]])
  std <- values(rst[["EFAamp"]])
  smx <- values(rst[["EFAsemax"]])
  
  #cor(cbind(avg,std,smx), method="spearman") %>% round(2)
  
  avgQt <- quantile(avg, qts, na.rm = TRUE)
  stdQt <- quantile(std, qts, na.rm = TRUE)
  
  recDF <- data.frame(avgCl = as.character(cut(avg, avgQt, labels = 1:4)),
                      stdCl = as.character(cut(std, stdQt, labels = 1:4)),
                      mmax  = as.character(smx))
  
  EFT_combns <- apply(recDF, MARGIN = 1, FUN = funConcat)
  
  message("\n\n-> Found ",length(unique(EFT_combns))," EFT combinations in year ",yr,"!\n\n")
  
  EFT_map <- raster(rst[[1]])
  values(EFT_map) <- EFT_combns
  
  #plot(EFT_map)
  
  fnOut <- paste(outFolder,"/LT578comp_EVI_jul-jun_EFT_",yr,"_med_amp_semax_025qts_v1.tif",sep="")
  writeRaster(EFT_map, fnOut, overwrite = TRUE)
  
}




