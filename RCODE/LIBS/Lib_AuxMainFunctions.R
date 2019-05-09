

## I/O Functions --------------------------------------------------------------------------



writeMultipleRasterFiles<-function(x,dates,outDir,prefix,suffix,...){
  
  if(!(is(x,"RasterStack") || is(x,"RasterBrick")))
    stop("Object in x must be of class RasterStack or RasterBrick")
  
  if(!dir.exists(outDir))
    stop("Output directory does not exists!")
  
  if(length(dates) != nlayers(x))
    stop("Number of layers in x must be equal to the number of dates ")
  
  ## Save raster files
  ##
  pb<-txtProgressBar(min = 1, max = nlayers(x), style = 3)
  
  for(i in 1:nlayers(x)){
    
    writeRaster(x[[i]],filename = paste(outDir,"/",prefix,".",dates[i],".",suffix,sep=""),...)
    
    setTxtProgressBar(pb, value=i)
  }
}




## Check objects class --------------------------------------------------------------------- 
# Raster time-series

isRasterBrickTS<-function(x) if(inherits(x,"RasterBrickTS")) return(TRUE) else return(FALSE)

isRasterStackTS<-function(x) if(inherits(x,"RasterStackTS")) return(TRUE) else return(FALSE)


## Dates processing in MODIS files --------------------------------------------------------

getDates<-function(x) gsub("A","",unlist(strsplit(x,"\\."))[2]) # Replaces the A prefix if exists


getMODISdates<-function(x) sapply(x, FUN = getDates)


MOD.getFiles<-function(x,ext="tif",...) list.files(x,pattern=paste(".",ext,"$",sep=""),...)


MOD.getYear<-function(x) as.integer(sapply(x,FUN=function(x) substr(gsub("A","",unlist(strsplit(x,"\\."))[2]),1,4)))


MOD.getDOY<-function(x) as.integer(sapply(x,FUN=function(x) substr(gsub("A","",unlist(strsplit(x,"\\."))[2]),5,7)))


# As YYYY-MM-DD
MOD.getDateChar<-function(x) sapply(x,FUN=function(x) as.character(as.Date(MOD.getDOY(x)-1,origin=paste(MOD.getYear(x),"01-01",sep="-"))))


MOD.getDate<-function(x) sapply(x,FUN=function(x) as.Date(MOD.getDOY(x)-1,origin=paste(MOD.getYear(x),"01-01",sep="-")))


getMODISfileDate<-function(x,which="all"){
  
  dt<-substr(unlist(strsplit(x,"\\."))[2],2,8)
  yr<-as.integer(substr(dt,1,4))
  jd<-as.integer(substr(dt,5,7))
  
  if(which=="all")
    return(c(dt,yr,jd))
  if(which=="dt")
    return(dt)
  if(which=="yr")
    return(yr)
  if(which=="jd")
    return(jd)
}


getYearSeq<-function(x){
  z<-unlist(strsplit(x,"/"))
  (as.integer(z[1])):(as.integer(z[2]))
} 



## Misc ------------------------------------------------------------------------------------------------------


randString <- function(n=12) paste(sample(c(sample(0:9,length(letters),replace=TRUE),sample(letters,length(letters))),n),collapse="")

# Creates a shapefile grid of raster centroids
# Uses expand.grid function and only works with projected CRS
#
getRasterCellXY<-function(x,outDir,shapeName,obj=FALSE,...){
  
  if(!dir.exists(outDir))
    stop("Output directory does not exists!")
  
  rst_resol<-res(x)[1]
  
  ext<-as.matrix(extent(x))
  x_min<-ext[1,1] 
  y_min<-ext[2,1]
  x_max<-ext[1,2] 
  y_max<-ext[2,2]
  
  x_seq<-seq(x_min+(rst_resol/2),x_max,250)
  y_seq<-seq(y_min+(rst_resol/2),y_max,250)
  
  xy_coord<-expand.grid(x_seq,y_seq)
  spPointsRstXY<-SpatialPointsDataFrame(xy_coord,proj4string = crs(x),data=data.frame(ID=1:nrow(xy_coord)))
  
  writeOGR(spPointsRstXY, outDir, shapeName, driver="ESRI Shapefile")
  
  if(obj){
    return(spPointsRstXY)
  }
  
}

# Creates a shapefile grid of raster centroids
# Uses rasterToPoints making it more general than getRasterCellXY

getGridFromRasterCentroids<-function(x,outDir,shapeName,obj=FALSE,...){
  
  if(!dir.exists(outDir))
    stop("Output directory does not exists!")
  
  if(is.character(x))
    x<-raster(x)

  pts<-rasterToPoints(x)
  spPtsDF<-SpatialPointsDataFrame(pts[,1:2],data=data.frame(ID=1:nrow(pts)),proj4string = crs(x))
  writeOGR(spPtsDF, outDir, shapeName, driver="ESRI Shapefile")
  
  if(obj){
    return(spPtsDF)
  }
  
}


## Calculate EF metrics --------------------------------------------------------------------------------------


# ''Productivity''
vi.mean<-function(x,...) mean(x,...)
vi.medv<-function(x,...) median(x,...)
vi.maxv<-function(x,...) max(x,...)
vi.minv<-function(x,...) min(x,...)
vi.qt05<-function(x,...) quantile(x,probs=0.05,...)
vi.qt25<-function(x,...) quantile(x,probs=0.25,...)
vi.qt75<-function(x,...) quantile(x,probs=0.75,...)
vi.qt95<-function(x,...) quantile(x,probs=0.95,...)


# ''Seasonality''
vi.stdv<-function(x,...) sd(x,...)
vi.covr<-function(x,...) sd(x,...)/mean(x,...)
vi.madv<-function(x,...) mad(x,...)

#vi.dqt1<-function(x,...) vi.qt95(x,...) - vi.qt05(x,...)
vi.dqt1<-function(x,...){ 
  qts<-quantile(x,probs=c(0.05,0.95))
  return(qts[2]-qts[1])
}

#vi.dqt2<-function(x,...) vi.qt75(x,...) - vi.qt25(x,...)

vi.dqt1n<-function(x,...) vi.qt95(x,...) - vi.qt05(x,...) / median(x,...)
vi.dqt2n<-function(x,...) vi.qt75(x,...) - vi.qt25(x,...) / median(x,...)


# ''Phenology/timming''
# vi.dmax<-function(x, ...){
#   which.max(x)
# }
# 
# vi.summ<-function(x,na.rm=FALSE){ # Convert to DOY and calculate the sin-transformation ("summerness")
#   doy<-seq(1,353,16)
#   sin((doy[x]/365)*pi) # Sine transform summerness
# } 
# 
# vi.wint<-function(x,na.rm=FALSE){
#   doy<-seq(1,353,16)
#   cos((doy[x]/365)*pi) # Cos transform summerness
# } 
# 

# This version of the function uses a time-series vector and calculates which (16-day) period
# is the maximum

# Springness/fall transform
vi.spring<-function(x,na.rm=FALSE){ 
  xi.max<-which.max(x)
  doy<-seq(1,353,16)
  sin((2*pi*(doy[xi.max]+11))/365)
} 

# Winter/summerness transform
vi.winter<-function(x,na.rm=FALSE){
  xi.max<-which.max(x)
  doy<-seq(1,353,16)
  cos((2*pi*(doy[xi.max]+11))/365)
} 





## Stability metrics ------------------------------------------------------------------------------------------


stlComponentsVar<-function(x,comp,FUN="CV",start,end,freq,...){
  
  x.ts<-ts(x,start=start,end=end,frequency=freq)
  stlDecomp<-stl(x.ts,...)
  
  if(FUN=="CV")
    return(sd(stlDecomp$time.series[,comp]) / mean(stlDecomp$time.series[,comp]))
  else if(FUN=="SD")
    return(sd(stlDecomp$time.series[,comp]))
  else if(FUN=="VAR")
    return(var(stlDecomp$time.series[,comp]))
  else
    stop("Invalid option in FUN!")
  
}

# 
# simpleHurstExp<- function(x) {
#   n <- length(x)
#   y <- x - mean(x)
#   s <- cumsum(y)
#   rs <- (max(s) - min(s))/sd(x)
#   log(rs)/log(n)
# }
# 
# 
# ## fArma: R/S Rescaled Range Statistic method
# 
# HurstExpRSfit<-function (x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5)){
#   
#   call = match.call()
#   data = list(x = x)
#   x = as.vector(x)
#   n = length(x)
#   increment = (log10(n/minnpts))/levels
#   M = floor(10^((1:levels) * increment))
#   M = M[M > 1]
#   Y = cumsum(x)
#   Y2 = cumsum(x * x)
#   RS = NULL
#   
#   for (m in M) {
#     S = sqrt(Y2[m]/m - (Y[m]/m)^2)
#     Z = Y[1:m] - (1:m) * Y[m]/m
#     STATS = (max(Z) - min(Z))/S
#     RS = c(RS, STATS)
# 
#   }
#   
#   wt = trunc((sign((M - cut.off[1]) * (cut.off[2] - M)) + 1)/2)
#   fit = lsfit(log10(M), log10(RS), wt)
#   return(fit$coef[[2]])
#   
# }
# 
# 
# 
#
# HurstExpRSfitCpp<-function (x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5)){
#   
#   call <- match.call()
#   data <- list(x = x)
#   x <- as.vector(x)
#   n <- length(x)
#   increment <- (log10(n/minnpts))/levels
#   M <- floor(10^((1:levels) * increment))
#   M <- M[M > 1]
#   Y <- cumsum(x)
#   Y2 <- cumsum(x * x)
#   #RS = NULL
#   # for (m in M) {
#   #   S = sqrt(Y2[m]/m - (Y[m]/m)^2)
#   #   Z = Y[1:m] - (1:m) * Y[m]/m
#   #   STATS = (max(Z) - min(Z))/S
#   #   RS = c(RS, STATS)
#   #   
#   # }
#   RS <- HurstFitRollingWindow(M, Y, Y2)
#   wt <- trunc((sign((M - cut.off[1]) * (cut.off[2] - M)) + 1)/2)
#   fit <- lsfit(log10(M), log10(RS), wt)
#   
#   return(fit$coef[[2]])
# }
# 
# 
# ## Corrected R over S Hurst exponent
# 
# HurstExpCorrectRoverS <- function(x) {
#   
#   # half intervals of indices
#   half <- function(N) sort(c(N, N[-length(N)]+((diff(N)+1)%/%2)))
#   
#   # define the R/S scale
#   rscalc <- function(x) {
#     n <- length(x)
#     y <- cumsum(x - mean(x))
#     R <- diff(range(y))
#     S <- sd(x)
#     return(R/S)
#   }
#   
#   # set initial values
#   X <- c(length(x))
#   Y <- c(rscalc(x))
#   N <- c(0, length(x) %/% 2, length(x))
#   
#   # compute averaged R/S for halved intervals
#   while ( min(diff(N)) >= 8 ) {
#     xl <- c()
#     yl <- c()
#     
#     for (i in 2:length(N)) {
#       rs <- rscalc(x[(N[i-1]+1):N[i]])
#       xl <- c(xl, N[i]-N[i-1])
#       yl <- c(yl, rs)
#     }
#     
#     X <- c(X, mean(xl))
#     Y <- c(Y, mean(yl))
#     
#     # next step
#     N <- half(N)
#   }
#   # apply linear regression
#   rs_lm <- lm(log(Y) ~ log(X))
#   
#   return(unname(coefficients(rs_lm)[2]))
# }
# 
# 

