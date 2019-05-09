

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

## ------------------------------------------------------------------------ ##


# Load Rcpp functions
#sourceCpp("D:/MyDocs/R-dev/PA_EcosystemStability/Cpp/HampelFilter.cpp")
#sourceCpp("D:/MyDocs/R-dev/PA_EcosystemStability/Cpp/setWeights.cpp")

# Load R auxiliary functions
#source("D:/MyDocs/R-dev/PA_EcosystemStability/R/LIBS/Lib_AuxMainFunctions.R")
#source("D:/MyDocs/R-dev/PA_EcosystemStability/R/LIBS/Lib_SmoothFunctions.R")



whitSmooth <- function(x, na.rm=TRUE, k=3, t.val=3.5, lambda = 4, w=NULL){
  
  # if(any(sapply(x, is.infinite))){
  #   return(NA)
  # }
  # 
  # x <- as.numeric(x)
  # #print(x)
  # 
  # if(na.rm & any(is.na(x))){
  #   return(NA)
  # }
  
  xHF <- try(hampel(x, k=k, t=t.val))
  
  if(inherits(xHF,"try-error")){
    return(rep(NA,length(x)))
  }else{
    if(is.list(xHF) & length(xHF)==2){
      x <- xHF[["y"]] # filtered values
    }else{
      return(rep(NA,length(x)))
    }
  }
  
  
  # x[x > 1] <- 1 # Replace spurious high-values
  # 
  # if(is.null(w))
  #   w <- rep(1,length(x))
  # 
  # sm <- suppressMessages(try(whit2(x, lambda = lambda, w = w), silent = TRUE))
  # 
  # if(any(sapply(sm, is.infinite))){
  #   return(NA)
  # }
  
  if(inherits(x,"try-error")){
    return(rep(NA,length(x)))
  }else{
    return(as.numeric(x))
  }
}




# Process data with raster calc and the pre-defined function ---------------------------------
#

fl <- list.files("./DATA/RASTER/S2_Vez/S2A_crop/NDVImax", full.names=TRUE)

fl <- fl[c(1:2,6:13,3:5)]

rstStack <- stack(fl)


#testPt1 <- matrix(c(544850, 4636060),1,2)
#testPt1 <- matrix(c(543208, 4635293),1,2) # Good!
#testPt1 <- matrix(c(546756, 4642453),1,2)
#testPt1 <- matrix(c(549331, 4637737),1,2)
#testPt1 <- matrix(c(546596, 4640707),1,2)
#testPt1 <- matrix(c(545049, 4640532),1,2)  # Floresta mista
#testPt1 <- matrix(c(544713, 4640881),1,2)  # Floresta carvalho? !Good!
#testPt1 <- matrix(c(543834, 4640760),1,2)  # Floresta mista?
#testPt1 <- matrix(c(543167, 4640311),1,2)  # Campo feno?
#testPt1 <- matrix(c(543174, 4640338),1,2)  # Campo feno?
#testPt1 <- matrix(c(558465, 4644045),1,2)  # Matos
#testPt1 <- matrix(c(556847, 4645211),1,2)  # Matos
testPt1 <- matrix(c(553674, 4646172),1,2)   # Floresta carvalho/mista?

pt1Data <- raster::extract(rstStack, testPt1) %>% as.numeric
#pt1Data[pt1Data > 2] <- 1
pt1Hf <- hampel(pt1Data, k=3, t=3.5)$y
pt1_Ws <- whitSmooth(pt1Data, lambda=4)
pt1_SGs <- savgol(pt1Hf, fl = 19, forder = 2)

plot(x=0:12, y=pt1Hf,type="p", xlab="Time index (month)", ylab="NDVI", ylim=c(0,1), pch=20, cex=2, 
     main = "NDVI values (Dec2015 - Dec2016)\nOak/mixed forest")
lines(x = 0:12, pt1_Ws, col="dark green", lwd=2.5)
lines(x = 0:12, pt1_SGs, lty=2, col="blue", lwd=2)
legend(x = "bottomright",legend=c("NDVI observed values","Whittaker smooth","Savitzky-Golay smooth"),
       col=c("black","dark green","blue"),
       lwd=c(NA, 2.5, 2),
       lty=c(NA, 1, 2),
       pch=c(19, NA, NA))


## ------------------------------------------------------------------------ ##

rstStack.whit2 <- raster::calc(rstStack,fun=whitSmooth)

