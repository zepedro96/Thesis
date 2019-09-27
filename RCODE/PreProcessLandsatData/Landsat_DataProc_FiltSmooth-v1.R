
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
library(sf)
library(dplyr)
library(ggnewscale)

# Load Rcpp functions
sourceCpp("./Cpp/setWeights.cpp")

# Load R auxiliary functions
source("./RCODE/LIBS/Lib_AuxMainFunctions.R")
source("./RCODE/LIBS/Lib_SmoothFunctions.R")




smoothTS_Hampel_Whittaker <- function(x, useHampelFilter=TRUE, HF.k=NULL, HF.t0=NULL, 
                                    lambda, useWt=TRUE, indices=NULL, quantProbs=NULL, quantWeights=NULL, 
                                    baseWeight=NULL, checkWts=TRUE, checkWin=3,na.rm=FALSE, impute.algo="ma",
                                    k.imp=4, st, en, fr){
  
  if(all(is.na(x))){
    return(NA)  
  }else{
    
    if(any(is.na(x))){
      x.ts <- ts(data = x, start = st, end = en, frequency = fr)
      suppressWarnings(x <- as.numeric(imputeTS::na.seadec(x.ts, algorithm="ma", k=k.imp)))
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
  indices      = c(rep(1,10),rep(2:35,each=12),rep(36,4))
  quantProbs   = c(0.5, 0.75, 0.9)
  quantWeights = c(2, 4, 10)
  baseWeight   = 1 
  checkWts     = FALSE
  checkWin     = 3
  na.rm        = na.rm
  st           = c(1984,3)
  en           = c(2019,4)
  fr           = 12
  
  return(smoothTS_Hampel_Whittaker(x, useHampelFilter, HF.k, HF.t0, #HF.winType, 
                                   lambda, useWt, indices, quantProbs, quantWeights, 
                                   baseWeight, checkWts, checkWin, na.rm ,st=st,en=en,fr=fr))
}

tsSmooth.HWwt.LT7 <- function(x, na.rm=TRUE){
  
  useHampelFilter = TRUE 
  HF.k  = 3 
  HF.t0 = 3.5 
  #HF.winType=1 
  lambda       = 5 
  useWt        = TRUE 
  indices      = c(rep(1:6,each=12), rep(7,4))
  quantProbs   = c(0.5, 0.75, 0.9)
  quantWeights = c(2, 4, 10)
  baseWeight   = 1 
  checkWts     = FALSE
  checkWin     = 3
  na.rm        = na.rm
  st           = c(2013,1)
  en           = c(2019,4)
  fr           = 12
  
  return(smoothTS_Hampel_Whittaker(x, useHampelFilter, HF.k, HF.t0, #HF.winType, 
                                   lambda, useWt, indices, quantProbs, quantWeights, 
                                   baseWeight, checkWts, checkWin, na.rm, st=st,en=en,fr=fr))
}


## ------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------ ##


LTdates <- readxl::read_excel("./DATA/RASTER/Landsat/EVI-32day/EVI-32-day_ImageList-v1.xlsx")
LT7_dates <- LTdates %>% filter(Satellite=="LE07", Year >= 2013)

LT7 <- brick("./DATA/RASTER/Landsat/EVI-32day/LT7_EVI_Vez_Basin_1999-2019.tif" ,values=TRUE)
LT7 <- LT7[[LT7_dates$ids]]
names(LT7) <- paste("EVIts_", LT7_dates$DateCode, sep="")


tsDates <- read.csv("./DATA/RASTER/Landsat/EVI-32day/FullImageList_1984-2019-LT-5-7-8_v2.csv")

fl <- list.files("./DATA/RASTER/Landsat/EVI-32day/full",pattern=".tif$", full.names = TRUE)
EVIts <- stack(fl)
names(EVIts) <- paste("EVIts_",tsDates$DateCode,sep="")


EVIts_whitsm <- calc(EVIts, fun = tsSmooth.HWwt)

# Write data to individual TIFF files
wdir <- "./DATA/RASTER/Landsat/EVI-32day/smoothed/"
N <- nlayers(EVIts_whitsm)
pb <- txtProgressBar(1,N,style = 3)

for(i in 1:N){
  
  dateCode <- DT_all$DateCode[i]
  satCode <- DT_all$Satellite[i]
  
  outFileName <- paste(wdir, "EVI_T1_32day_WhitSmooth_",dateCode,"_",satCode,"-v1.tif",sep="")
  
  if(!file.exists(outFileName)){
    cat("\n-> Writing file:",outFileName,"\n\n")
    writeRaster(EVIts_whitsm[[i]],outFileName)
  }else{
    cat("\n-> Skipping file!!\n\n")
  }
  setTxtProgressBar(pb, i)
}


## ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- ##


#testPt <- matrix(c(544713, 4640881),1,2)  # Floresta carvalho? !Good!
testPt <- matrix(c(553674, 4646172),1,2)   # Floresta carvalho/mista?

pt1Data <- raster::extract(EVIts, testPt) %>% as.numeric
pt1Data_LT7 <- raster::extract(LT7, testPt) %>% as.numeric

ptdataTs <- ts(data = pt1Data, start = c(1984,3), end=c(2019,4), frequency = 12)
ptdataTs_LT7 <- ts(data = pt1Data_LT7, start = c(2013,1), end=c(2019,4), frequency = 12)

time(ptdataTs)
time(ptdataTs_LT7)

plot(ptdataTs, type = "l")
plot(ptdataTs_LT7, type = "l")

pt1DataImpute_ma <- imputeTS::na.seadec(ptdataTs, algorithm = "ma", k = 4)
plot(pt1DataImpute_ma, type = "l")

pt1DataImpute_ma_LT7 <- imputeTS::na.seadec(ptdataTs_LT7, algorithm = "ma", k = 4)
plot(pt1DataImpute_ma_LT7, type = "l")


smWhit     <- tsSmooth.HWwt(pt1Data)
smWhit_LT7 <- tsSmooth.HWwt.LT7(pt1Data_LT7)


LT_DF <- data.frame( EVI    = pt1Data, 
                     EVI_   = pt1DataImpute_ma, 
                     smType = "Whittaker",
                     EVI_sm = smWhit, 
                     dates  = as.Date(tsDates$ImgDate), 
                     sat    = tsDates$Satellite)



LT7_DF <- data.frame(EVI    = pt1Data_LT7, 
                     EVI_   = pt1DataImpute_ma_LT7, 
                     smType = "Whittaker",
                     EVI_sm = smWhit_LT7, 
                     dates  = as.Date(LT7_dates$ImgDate), 
                     sat    = LT7_dates$Satellite)


vals <-  c("Composite series LT 5,7,8 (sm.)" = "red",
          "Landsat 7 ETM+ (sm.)" = "black")
vals1 <- c("LT05" = "blue",
          "LE07" = "orange",
          "LC08" = "dark green")


g <- ggplot(LT_DF)+
  geom_line(aes(x = dates, y = EVI_),linetype="dashed", size=0.5, color="black") +
  geom_point(aes(x = dates, y = EVI_, shape = sat, color=sat), size=2) +
   
  #new_scale_color() +
  #scale_colour_manual(name = "sat", values=vals1 )
  geom_line(aes(x = dates, y = EVI_sm, color="Composite series LT 5,7,8 (sm.)"), size=1.25) + 
  geom_line(data=LT7_DF, mapping=aes(x = dates, y = EVI_sm, color="LE07"), size=1.25) +
  #scale_colour_manual(name = "Smoothed series:", values=vals)
  scale_color_brewer(palette = "Set1", name="Data source:") +
  scale_shape(name="") +
  theme_bw() + 
  theme(legend.position="bottom") + 
  ylab("EVI (Enhanced Vegetation Index)") + 
  xlab("Time/Year")

plot(g)


## ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- ##


## Read base data
##
LTdates <- readxl::read_excel("./DATA/RASTER/Landsat/EVI-32day/EVI-32-day_ImageList-v1.xlsx")
LT7_dates <- LTdates %>% filter(Satellite=="LE07", Year >= 2013)

LT7 <- brick("./DATA/RASTER/Landsat/EVI-32day/LT7_EVI_Vez_Basin_1999-2019.tif" ,values=TRUE)
LT7 <- LT7[[LT7_dates$ids]]
names(LT7) <- paste("EVIts_", LT7_dates$DateCode, sep="")


tsDates <- read.csv("./DATA/RASTER/Landsat/EVI-32day/FullImageList_1984-2019-LT-5-7-8_v2.csv")

fl <- list.files("./DATA/RASTER/Landsat/EVI-32day/full",pattern=".tif$", full.names = TRUE)
EVIts <- stack(fl)
names(EVIts) <- paste("EVIts_",tsDates$DateCode,sep="")



testPoints <- read_sf("./DATAtoShare/testPoints.shp") %>% select(Name, Name_EN) %>% st_zm %>% 
  st_transform(crs=c("+proj=utm +zone=29 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# ard <- read_sf("./DATAtoShare/area ardida.shp") %>% select(Name) %>% st_zm %>%
#   st_transform(crs=c("+proj=utm +zone=29 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

#pts <- (rbind(testPoints,ard))
pts_DF <- `st_geometry<-`(testPoints, NULL)



for(i in 26:nrow(pts)){
  
  
  ptType <- pts_DF[i,"Name_EN",drop=TRUE]
  
  pt1Data     <- raster::extract(EVIts, pts[i,]) %>% as.numeric
  pt1Data_LT7 <- raster::extract(LT7, pts[i,]) %>% as.numeric
  
  ptdataTs             <- ts(data = pt1Data, start = c(1984,3), end=c(2019,4), frequency = 12)
  ptdataTs_LT7         <- ts(data = pt1Data_LT7, start = c(2013,1), end=c(2019,4), frequency = 12)
  pt1DataImpute_ma     <- imputeTS::na_seadec(ptdataTs, algorithm = "ma", k = 4)
  pt1DataImpute_ma_LT7 <- imputeTS::na_seadec(ptdataTs_LT7, algorithm = "ma", k = 4)
  
  
  smWhit     <- tsSmooth.HWwt(pt1Data)
  smWhit_LT7 <- tsSmooth.HWwt.LT7(pt1Data_LT7)
  
  
  LT_DF <- data.frame( EVI    = pt1Data, 
                       EVI_   = pt1DataImpute_ma, 
                       smType = "Whittaker",
                       EVI_sm = smWhit, 
                       dates  = as.Date(tsDates$ImgDate), 
                       sat    = tsDates$Satellite)
  
  LT7_DF <- data.frame(EVI    = pt1Data_LT7, 
                       EVI_   = pt1DataImpute_ma_LT7, 
                       smType = "Whittaker",
                       EVI_sm = smWhit_LT7, 
                       dates  = as.Date(LT7_dates$ImgDate), 
                       sat    = LT7_dates$Satellite)
  
  g <- ggplot(LT_DF)+
    geom_line(aes(x = dates, y = EVI_),linetype="dashed", size=0.5, color="black") +
    geom_point(aes(x = dates, y = EVI_, shape = sat, color=sat), size=2.5) +
    geom_line(aes(x = dates, y = EVI_sm, color="Composite series LT 5,7,8\n(Whittaker smooth)"), size=1.25) + 
    #geom_line(data=LT7_DF, mapping=aes(x = dates, y = EVI_sm, color="LE07"), size=1.25) +
    scale_color_brewer(palette = "Set1", name="Data source:") +
    scale_shape(name="") +
    ylab("EVI (Enhanced Vegetation Index)") + 
    xlab("Time/Year") + 
    theme_bw() + 
    theme(legend.position="bottom") + 
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.title = element_text(size=22),
          plot.subtitle = element_text(size=20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=20)) +
    labs(title = "EVI 32-day composite time series for Landsat 5,7,8 (1984 - 2019)",
         subtitle = paste("LULC type: ",ptType," [",i,"]",sep=""))
  
  plot(g)
  
  ggplot2::ggsave(plot=g, filename=paste("./OUTtoShare/tsPlots_v2/EVI_32day_LTcompositeTS_[",i,"]_",
                                         gsub("\\/","_",ptType),".png",sep=""),
                  width = 14, height = 8)
  
  cat("\n---> Finished point",i,"out of",nrow(pts),"\n\n")
  
}




