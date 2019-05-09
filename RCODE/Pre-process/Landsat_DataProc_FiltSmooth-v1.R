
library(Rcpp)
library(sp)
library(raster)
library(rts)
library(rgdal)
library(SparseM)
library(ptw)
library(pracma)
library(magrittr)
library(ggplot2)


# Load Rcpp functions
sourceCpp("./Cpp/setWeights.cpp")

# Load R auxiliary functions
source("./RCODE/LIBS/Lib_AuxMainFunctions.R")
source("./RCODE/LIBS/Lib_SmoothFunctions.R")




smoothTS_Hampel_Whittaker <- function(x, useHampelFilter=TRUE, HF.k=NULL, HF.t0=NULL, 
                                    lambda, useWt=TRUE, indices=NULL, quantProbs=NULL, quantWeights=NULL, 
                                    baseWeight=NULL, checkWts=TRUE, checkWin=3,na.rm=FALSE, impute.algo="ma",
                                    k.imp=4){
  
  if(all(is.na(x))){
    return(NA)  
  }else{
    
    if(any(is.na(x))){
      suppressWarnings(x <- imputeTS::na.seadec(x, algorithm="ma", k=k.imp))
    }
    
    ## Use Hampel filter/identifier
    if(useHampelFilter){
      if(is.null(HF.k) || is.null(HF.t0)){
        stop("Need to set k and t0 for use in Hampel filter/identifier!")
      }
    }else{
      #x<-HampelFilterCpp(x, HF.k, HF.t0, HF.winType)
      x <- hampel(x, HF.k, HF.t0)$y
    } 
    ## Use weights
    if(useWt){
      if(is.null(indices) || is.null(quantProbs) || is.null(quantWeights) || is.null(baseWeight))
        stop("Parameters indices, quantProbs, quantWeights or baseWeight can not be set to NULL if useWt is TRUE!")
      ## Rcpp function to weight
      w <- setWeightsCpp(x,indices, quantProbs, quantWeights, baseWeight, checkWts, checkWin) 
      return(whit2(x, lambda, w))
    }else{
      return(whit2(x,lambda))
    }
  }
}

tsSmooth.HWwt <- function(x, na.rm=TRUE){
  
  useHampelFilter = TRUE 
  HF.k  = 3 
  HF.t0 = 3.5 
  #HF.winType=1 
  lambda       = 5 
  useWt        = TRUE 
  indices      = c(rep(1,9),rep(2:6,each=12),rep(7,4))
  quantProbs   = c(0.5, 0.75, 0.9)
  quantWeights = c(2, 5, 10)
  baseWeight   = 1 
  checkWts     = FALSE
  checkWin     = 3
  na.rm        = na.rm
  
  return(smoothTS_Hampel_Whittaker(x, useHampelFilter, HF.k, HF.t0, #HF.winType, 
                                   lambda, useWt, indices, quantProbs, quantWeights, 
                                   baseWeight, checkWts, checkWin, na.rm))
  
}


NDVIts <- stack("./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_Vez_Basin_2013-2019.tif")


NDVIts_smWhitt <- calc(NDVIts, fun = tsSmooth.HWwt)

writeRaster(NDVIts_smWhitt,filename = "./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_Vez_Basin_2013_2019_smWhitt.tif")


ptDataOriginal <- raster::extract(NDVIts, testPt) %>% as.numeric
ptDataSmooth <- raster::extract(NDVIts_smWhitt, testPt) %>% as.numeric




# ## ------------------------------------------------------------------------ ##
# 
# tsDates <- #rbind(
#   #readxl::read_excel("./DATA/RASTER/Landsat/NDVI-32day/Landsat_5-8_NDVI-32Day_imageDates.xlsx", 1),
#   readxl::read_excel("./DATA/RASTER/Landsat/NDVI-32day/Landsat_5-8_NDVI-32Day_imageDates.xlsx", 2)#)
# 
# NDVIts <- stack(#c("./DATA/RASTER/Landsat/NDVI-32day/LT5_NDVI_Vez_Basin_1990-2012.tif",
#                   "./DATA/RASTER/Landsat/NDVI-32day/LT8_NDVI_Vez_Basin_2013-2019.tif")
# names(NDVIts) <- paste("NDVI_",tsDates$DateCode,sep="")
# 
# 
# #testPt <- matrix(c(544713, 4640881),1,2)  # Floresta carvalho? !Good!
# testPt <- matrix(c(553674, 4646172),1,2)   # Floresta carvalho/mista?
# pt1Data <- raster::extract(NDVIts, testPt) %>% as.numeric
# 
# #pt1DataSmWhitt <- raster::extract(NDVIts_smWhitt, testPt) %>% as.numeric
# 
# pt1DataImpute <- imputeTS::na.seadec(pt1Data, algorithm="ma", k=4)
# 
# # qt <- quantile(pt1DataImpute, 0.75)
# # w <- rep(1, length(pt1Data))
# # w[pt1DataImpute > qt] <- 5
# # smWhit <- whit2(pt1DataImpute, lambda=4, w = w)
# 
# yrs <- tsDates$Year
# names(yrs) <- tsDates$Year
# yr <- sort(unique(yrs))
# yr_ <- 1:length(yr)
# names(yr_) <- yr
# inds <- yr_[yrs]
# 
# 
# 
# smWhit <- 
#   smoothTS_Hampel_Whittaker(pt1DataImpute, useHampelFilter=TRUE, HF.k=3, HF.t0=3, 
#                             lambda=5, useWt=TRUE, indices=inds, quantProbs=c(0.05, 0.5, 0.75, 0.9), 
#                             quantWeights=c(0.1, 1, 2, 10), 
#                             baseWeight=1, checkWts=FALSE, checkWin=3,na.rm=FALSE)
# 
# smSavGol <- savgol(pt1DataImpute, fl = 9, forder = 2)
# 
# testDF <- rbind(data.frame(ndvi=pt1Data, ndvi_ = pt1DataImpute, smType = "Whittaker", 
#                            ndvi_sm = smWhit, dates=as.Date(tsDates$ImgDate))
#                 # 
#                 # data.frame(ndvi=pt1Data, ndvi_ = pt1DataImpute, smType = "Savitzky-Golay", 
#                 #            ndvi_sm = smSavGol, dates=as.Date(tsDates$ImgDate))
#                 
#                 )
# 
# g <- ggplot(testDF)+
#   geom_line(aes(x = dates, y=ndvi_),linetype="dashed") + 
#   geom_point(aes(x = dates, y=ndvi_), size=3, color="red") + 
#   geom_line(aes(x = dates, y = ndvi_sm, color=smType), size=1.25) 
#   # scale_color_manual(name="Smoothing\nalgorithm:",
#   #                    values=c("Savitzky-Golay"="blue", "Whittaker"="dark green"))
# 
# plot(g)
# 
# 
