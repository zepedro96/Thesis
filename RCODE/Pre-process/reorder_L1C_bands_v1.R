
library(raster)
library(utils)
library(tools)

list.files("F:/DATA/L1C/S2A_MSIL1C_20160213T113202_N0201_R080_T29TNG_20160213T113807.SAFE")

fl <- list.files("D:/MyDocs/R-dev/Thesis/DATA/RASTER/S2_Vez/S2A_crop/SpectralBands", pattern = ".SAFE_",
                 full.names=TRUE, recursive=FALSE)


outFolder <- "D:/MyDocs/R-dev/Thesis/DATA/RASTER/S2_Vez/S2A_crop/SpectralBands"


for(i in 1:length(fl)){
  
  r <- stack(fl[i])
  fnOut <- paste(file_path_sans_ext(basename(fl[i])),"_v2-reord.tif",sep="")
  
  # Bands order: b2, b3, b4, b8 (10m) / b11, b12, b5, b6, b7, b8a (20m)
  r1 <- r[[c(1,2,3,7,8,9,4,5,6,10)]]
 
  writeRaster(r1, filename = paste(outFolder,"/",fnOut,sep="")) 
  print(i)
}

