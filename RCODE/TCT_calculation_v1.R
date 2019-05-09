
library(raster)
library(tools)
library(utils)

# Bands list: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b8a
#
# brCoeffs <- c(0.0356, 0.0822, 0.1360, 0.2611, 0.2964, 0.3338, 0.3877, 
#              0.3895, 0.949, 0.0009, 0.3882, 0.1366, 0.4750)
# 
# grCoeffs <- c(-0.0635, -0.1128, -0.1680, -0.3480, -0.3303, 0.0852, 
#              0.3302, 0.3165, 0.0467, -0.0009, -0.4578, -0.4064, 0.3625)
# 
# wtCoeffs <- c(0.0649, 0.1363, 0.2802, 0.3072, 0.5288, 0.1379, -0.0001,
#               -0.0807, -0.0302, 0.0003, -0.4064, -0.5602, -0.1389)
# 


#bandIndices <- c(1,2,3,4,5,6) # --> raster to band order: b2, b3, b4, b8, b11 and b12

# b2, b3, b4, b8, b11, b12
wt <- c(0.1509, 0.1973, 0.3279, 0.3406, -0.7112, 0.4572)

# b2, b3, b4, b8, b10, b12 (replaced b10 by-> b11)
br <- c(0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.186312)

# b2, b3, b4, b8, b11, b12
gr <- c(-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800)


fl <- list.files("D:/MyDocs/R-dev/Thesis/DATA/RASTER/S2_Vez/S2A_crop/SpectralBands",
                 full.names=TRUE, recursive=FALSE)

pb <- txtProgressBar(1,length(fl),style=3)

for(i in 1:length(fl)){
  
  fnOut <- paste(file_path_sans_ext(basename(fl[i])),"_TCTwgb.tif",sep="")
  r <- (stack(fl[i])) * 1E-4 
    
  wtRst <- (wt[1]*r[[1]]) + (wt[2]*r[[2]]) + (wt[3]*r[[3]]) + (wt[4]*r[[4]]) + (wt[5]*r[[5]]) + (wt[6]*r[[6]])
  brRst <- (br[1]*r[[1]]) + (br[2]*r[[2]]) + (br[3]*r[[3]]) + (br[4]*r[[4]]) + (br[5]*r[[5]]) + (br[6]*r[[6]])
  grRst <- (gr[1]*r[[1]]) + (gr[2]*r[[2]]) + (gr[3]*r[[3]]) + (gr[4]*r[[4]]) + (gr[5]*r[[5]]) + (gr[6]*r[[6]])
  
  #plot(wtRst)
   
  writeRaster(stack(wtRst, grRst, brRst), 
              filename=paste("./DATA/RASTER/S2_Vez/S2A_crop/TCT/",fnOut,sep=""))
  
  setTxtProgressBar(pb,i)
  
}

