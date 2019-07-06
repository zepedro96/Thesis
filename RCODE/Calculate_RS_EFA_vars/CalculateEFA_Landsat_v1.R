
library(raster)
library(dplyr)
library(readxl)

tsDates <- read.csv("./DATA/RASTER/Landsat/EVI-32day/FullImageList_1984-2019-LT-5-7-8_v2.csv")

fl <- list.files("./DATA/RASTER/Landsat/EVI-32day/smoothed", pattern=".tif$", full.names = TRUE)
EVIts_sm <- stack(fl)
names(EVIts_sm) <- paste("EVIts_",tsDates$DateCode,sep="")

## ----------------------------------------------------------------------------------------------------- ##


vi.dqt1<-function(x,...){ 
  qts<-quantile(x,probs=c(0.1,0.90))
  return(qts[2]-qts[1])
}

vi.spring<-function(x,na.rm=FALSE){ 
  xi.max<-which.max(x)
  sin((2*pi*xi.max)/12)
} 

# Winter/summerness transform
vi.winter<-function(x,na.rm=FALSE){
  xi.max<-which.max(x)
  cos((2*pi*xi.max)/12)
} 

vi.monthMax <- function(x,na.rm=FALSE){
  return(which.max(x))
}

vi.monthMax2 <- function(x, na.rm=FALSE){
  if(length(x)==4){ # ts start 1984
    mts <- 3:6
  }else if(length(x)==10){
    mts <- c(7:12,1:4) # ts end 2019
  }else{
    mts <- c(7:12,1:6)
  }
  return(mts[which.max(x)])
}

## ----------------------------------------------------------------------------------------------------- ##

# January - December seq ids
# inds <- c(rep(1,10),rep(2:35,each=12),rep(36,4))
# length(inds)

# July - June seq ids
inds <- c(rep(1,4),rep(2:36,each=12))[1:nlayers(EVIts_sm)]
length(inds)

monthsSeq <- c(3:12,rep(1:12,36))[1:nlayers(EVIts_sm)]
indsDF <- as.data.frame(cbind(index = inds,
                              Yr    = tsDates$Year, 
                              Month = tsDates$Month))

print(indsDF)

indsDF %>% group_by(index) %>% summarize(n=n())

## ----------------------------------------------------------------------------------------------------- ##

EFA_min <- stackApply(EVIts_sm, indices = inds, fun = min)
EFA_max <- stackApply(EVIts_sm, indices = inds, fun = max)
EFA_dqt <- stackApply(EVIts_sm, indices = inds, fun = vi.dqt1)
EFA_amp <- EFA_max - EFA_min
EFA_med <- stackApply(EVIts_sm, indices = inds, fun = median)

EFA_avg <- stackApply(EVIts_sm, indices = inds, fun = mean)
EFA_std <- stackApply(EVIts_sm, indices = inds, fun = sd)
EFA_cfv <- EFA_std / EFA_avg
ESPI <- EFA_avg - EFA_std

EFA_sprg <- stackApply(EVIts_sm, indices = inds, fun = vi.spring)
EFA_wint <- stackApply(EVIts_sm, indices = inds, fun = vi.winter)


EFA_mmax <- stackApply(EVIts_sm, indices = inds, fun = vi.monthMax2)

reclMat <- cbind(is=1:12, becomes=c(1,1,2,2,2,3,3,3,4,4,4,1))
EFA_mmax_seas <- reclassify(EFA_mmax, reclMat)

plot(EFA_mmax_seas[[31]])

## ----------------------------------------------------------------------------------------------------- ##


writeRaster(EFA_min,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAmin.tif")
writeRaster(EFA_max,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAmax.tif")
writeRaster(EFA_med,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAmed.tif")

writeRaster(EFA_amp,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAamp.tif")
writeRaster(EFA_dqt,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAdqt.tif")

writeRaster(EFA_sprg,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAsprg.tif")
writeRaster(EFA_wint,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAwint.tif")

writeRaster(EFA_mmax,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAmmax.tif", overwrite=TRUE)
writeRaster(EFA_mmax_seas,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAsemax.tif", overwrite=TRUE)

writeRaster(EFA_avg,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAavg.tif")
writeRaster(EFA_std,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAstd.tif")
writeRaster(EFA_cfv,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAcv.tif")
writeRaster(ESPI,"./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_ESPI.tif")

## ----------------------------------------------------------------------------------------------------- ##


EFA_min  <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAmin.tif")
EFA_max  <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAmax.tif")
EFA_med  <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAmed.tif")

EFA_amp  <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAamp.tif")
EFA_dqt  <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAdqt.tif")

EFA_sprg <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAsprg.tif")
EFA_wint <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAwint.tif")

EFA_mmax <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAmmax.tif")
EFA_mmax_seas <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAsemax.tif")

EFA_avg <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAavg.tif")
EFA_std <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAstd.tif")
EFA_cfv <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_EFAcv.tif")
ESPI    <- stack("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/LT578comp_EVI_1984_2019_jul-jun_ESPI.tif")










