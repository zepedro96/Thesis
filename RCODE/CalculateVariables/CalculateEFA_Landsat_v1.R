

vi.dqt1<-function(x,...){ 
  qts<-quantile(x,probs=c(0.05,0.95))
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


inds <- c(rep(1,9),rep(2:6,each=12),rep(7,4))

EFA_min <- stackApply(NDVIts_smWhitt, indices = inds, fun = min)
EFA_max <- stackApply(NDVIts_smWhitt, indices = inds, fun = max)
EFA_amp <- EFA_max - EFA_min
EFA_med <- stackApply(NDVIts_smWhitt, indices = inds, fun = median)

EFA_avg <- stackApply(NDVIts_smWhitt, indices = inds, fun = mean)
EFA_std <- stackApply(NDVIts_smWhitt, indices = inds, fun = sd)
EFA_cfv <- EFA_std / EFA_avg
ESPI <- EFA_avg - EFA_std

EFA_sprg <- stackApply(NDVIts_smWhitt, indices = inds, fun = vi.spring)
EFA_wint <- stackApply(NDVIts_smWhitt, indices = inds, fun = vi.winter)


writeRaster(EFA_min,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_EFAmin.tif")
writeRaster(EFA_max,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_EFAmax.tif")
writeRaster(EFA_med,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_EFAmed.tif")
writeRaster(EFA_amp,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_EFAamp.tif")
writeRaster(EFA_sprg,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_EFAsprg.tif")
writeRaster(EFA_wint,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_EFAwint.tif")

writeRaster(EFA_avg,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_EFAavg.tif")
writeRaster(EFA_std,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_EFAstd.tif")
writeRaster(EFA_cfv,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_EFAcv.tif")
writeRaster(ESPI,"./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_2013-19_ESPI.tif")