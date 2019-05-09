


## João Gonçalves
## Porto, 2016
##
## Ecosystem stability analyses of Portugal protected areas
##
## Simplifies several functions from the earlywarnings package to 
## allow their application over raster datasets
##

library(earlywarnings)


generic_ews_mod <- function(timeseries, winsize = 50, detrending = c("no", "gaussian", 
    "loess", "linear", "first-diff"), bandwidth = NULL, span = NULL, degree = NULL, 
    logtransform = FALSE, interpolate = FALSE, AR_n = FALSE, powerspectrum = FALSE, 
    ar1=TRUE, sd=TRUE, sk=TRUE, kurt=TRUE, cv=TRUE, returnrate=TRUE, densratio=TRUE, 
    acf1=TRUE, outputByTimeIndex=TRUE) {
    
    # timeseries<-ts(timeseries)
    timeseries <- as.matrix(timeseries)  #strict data-types the input data as tseries object for use in later steps
    if (dim(timeseries)[2] == 1) {
        Y = timeseries
        timeindex = 1:dim(timeseries)[1]
    } else if (dim(timeseries)[2] == 2) {
        Y <- timeseries[, 2]
        timeindex <- timeseries[, 1]
    } else {
        warning("not right format of timeseries input")
    }
    
    # Interpolation
    if (interpolate) {
        YY <- approx(timeindex, Y, n = length(Y), method = "linear")
        Y <- YY$y
    } else {
        Y <- Y
    }
    
    # Log-transformation
    if (logtransform) {
        Y <- log(Y + 1)
    }
    
    # Detrending
    detrending <- match.arg(detrending)
    if (detrending == "gaussian") {
        if (is.null(bandwidth)) {
            bw <- round(bw.nrd0(timeindex))
        } else {
            bw <- round(length(Y) * bandwidth/100)
        }
        smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = bw, range.x = range(timeindex), 
            x.points = timeindex)
        nsmY <- Y - smYY$y
        smY <- smYY$y
    } else if (detrending == "linear") {
        nsmY <- resid(lm(Y ~ timeindex))
        smY <- fitted(lm(Y ~ timeindex))
    } else if (detrending == "loess") {
        if (is.null(span)) {
            span <- 25/100
        } else {
            span <- span/100
        }
        if (is.null(degree)) {
            degree <- 2
        } else {
            degree <- degree
        }
        smYY <- loess(Y ~ timeindex, span = span, degree = degree, normalize = FALSE, 
            family = "gaussian")
        smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
        nsmY <- Y - smY
    } else if (detrending == "first-diff") {
        nsmY <- diff(Y)
        timeindexdiff <- timeindex[1:(length(timeindex) - 1)]
    } else if (detrending == "no") {
        smY <- Y
        nsmY <- Y
    }
    
    # Rearrange data for indicator calculation
    mw <- round(length(Y) * winsize/100)
    omw <- length(nsmY) - mw + 1  ##number of moving windows
    low <- 2
    high <- mw
    nMR <- matrix(data = NA, nrow = mw, ncol = omw)
    x1 <- 1:mw
    
    for (i in 1:omw) {
        Ytw <- nsmY[i:(i + mw - 1)]
        nMR[, i] <- Ytw
    }

    # Calculate indicators
    nARR <- numeric()
    nSD <- numeric()
    nSK <- numeric()
    nKURT <- numeric()
    nACF <- numeric()
    nDENSITYRATIO <- numeric()
    nSPECT <- matrix(0, nrow = mw, ncol = ncol(nMR))
    nCV <- numeric()
    smARall <- numeric()
    smARmaxeig <- numeric()
    detB <- numeric()
    ARn <- numeric()
    
    if(sd) nSD <- apply(nMR, 2, sd, na.rm = TRUE)
    nn<-ncol(nMR)
    
    for (i in 1:nn) {
      if(ar1) nYR <- ar.ols(nMR[, i], aic = FALSE, order.max = 1, demean = TRUE, intercept = FALSE)
      if(ar1) nARR[i] <- nYR$ar
      if(sk) nSK[i] <- abs(moments::skewness(nMR[, i], na.rm = TRUE))
      if(kurt) nKURT[i] <- moments::kurtosis(nMR[, i], na.rm = TRUE)
      if(cv) nCV[i] <- nSD[i]/mean(nMR[, i])
      if(acf1) ACF <- acf(nMR[, i], lag.max = 1, type = c("correlation"), plot = FALSE)
      if(acf1) nACF[i] <- ACF$acf[2]
      if(densratio) spectfft <- spec.ar(nMR[, i], n.freq = mw, plot = FALSE, order = 1)
      if(densratio) nSPECT[, i] <- spectfft$spec
      if(densratio) nDENSITYRATIO[i] <- spectfft$spec[low]/spectfft$spec[high]
      
      if (AR_n) {
        ## RESILIENCE IVES 2003 Indicators based on AR(n)
        ARall <- ar.ols(nMR[, i], aic = TRUE, order.max = 6, demean = F, intercept = F)
        smARall[i] <- ARall$ar[1]
        ARn[i] <- ARall$order
        roots <- Mod(polyroot(c(rev(-ARall$ar), 1)))
        smARmaxeig[i] <- max(roots)
        detB[i] <- (prod(roots))^(2/ARn[i])
      }
    }
    
    if(returnrate) nRETURNRATE = 1/nARR
    
    # Estimate Kendall trend statistic for indicators
    timevec <- seq(1, nn)
    
    if(ar1) KtAR <- cor.test(timevec, nARR, alternative = c("two.sided"), method = c("kendall"))$estimate
    if(acf1) KtACF <- cor.test(timevec, nACF, alternative = c("two.sided"), method = c("kendall"))$estimate
    if(sd) KtSD <- cor.test(timevec, nSD, alternative = c("two.sided"), method = c("kendall"))$estimate
    if(sk) KtSK <- cor.test(timevec, nSK, alternative = c("two.sided"), method = c("kendall"))$estimate
    if(kurt) KtKU <- cor.test(timevec, nKURT, alternative = c("two.sided"), method = c("kendall"))$estimate
    if(densratio) KtDENSITYRATIO <- cor.test(timevec, nDENSITYRATIO, alternative = c("two.sided"))$estimate
    if(returnrate) KtRETURNRATE <- cor.test(timevec, nRETURNRATE, alternative = c("two.sided"))$estimate
    if(cv) KtCV <- cor.test(timevec, nCV, alternative = c("two.sided"), method = c("kendall"))$estimate
    
    
    outDF <- data.frame(timeindex[mw:length(nsmY)], row.names = NULL)
    cn <- c("timeindex")
    kts<-c()
    
    if(ar1){
      outDF <- data.frame(outDF,nARR)
      cn <- c(cn, "ar1")  
      kts<-c(kts, KtAR)
    }
    
    if(sd){
      outDF <- data.frame(outDF,nSD)
      cn <- c(cn, "sd")  
      kts<-c(kts, KtSD)
    }
    
    if(sk){
      outDF <- data.frame(outDF,nSK)
      cn <- c(cn, "sk")  
      kts<-c(kts, KtSK)
    }
    
    if(kurt){
      outDF <- data.frame(outDF,nKURT)
      cn <- c(cn, "kurt")  
      kts<-c(kts, KtKU)
    }
    
    if(cv){
      outDF <- data.frame(outDF,nCV)
      cn <- c(cn, "cv")  
      kts<-c(kts, KtCV)
    }
    
    if(returnrate){
      outDF <- data.frame(outDF,nRETURNRATE)
      cn <- c(cn, "returnrate")  
      kts<-c(kts, KtRETURNRATE)
    }
    
    if(densratio){
      outDF <- data.frame(outDF,nDENSITYRATIO)
      cn <- c(cn, "densratio")  
      kts<-c(kts, KtDENSITYRATIO)
    }
    
    if(acf1){
      outDF <- data.frame(outDF,nACF)
      cn <- c(cn, "acf1")  
      kts<-c(kts, KtACF)
    }
    
    colnames(outDF) <- cn
    names(kts) <- cn[-1]
    
    # Output
    if(outputByTimeIndex){
      return(list(ews_metrics = outDF, KendallTau = kts))
    }else{
      return(kts)
    }
} 


ews_rw<-function(x, winsize = 20, na.rm = TRUE){
  
  # Use a mask in the first layer/or x value
  if(x[1]==0)
    return(NA)
  else{
    # Make ts object
    x.ts<-ts(data=x[-1],start=c(2000,4),end=c(2016,16),frequency=23)
    
    res<-try(generic_ews_mod(x.ts,winsize,detrending = "no", 
                             bandwidth = NULL, span = NULL, degree = NULL, 
                             logtransform = FALSE, interpolate = FALSE, 
                             AR_n = FALSE, powerspectrum = FALSE, ar1=TRUE, 
                             sd=TRUE, sk=TRUE, kurt=TRUE, cv=TRUE, returnrate=TRUE, 
                             densratio=FALSE, acf1=TRUE, outputByTimeIndex=FALSE), 
             silent = TRUE)
    
    if(inherits(res,"try-error"))
      return(NA)
    else
      return(res)
    }

}




calculate_AR_mw<-function(Yt,mw_size=30,mw_type="perc",type="ACF",out="trend_kendall",...){
  
  Yt<-ts(data=Yt,start=c(2000,4),end=c(2016,16),frequency=23)
  
  # if(!is.ts(Yt))
  #   stop("Yt must be a ts object")
  
  # Rearrange data for indicator calculation
  if(mw_type=="perc") mw <- round(length(Yt) * mw_size/100)
  if(mw_type=="npts") mw <- mw_size
  
  omw <- length(Yt) - mw + 1  ##number of moving windows
  ewsOut<-vector(mode = "numeric", length = omw)

  
  for (i in 1:omw) {
    
    # data in window
    Ytw <- Yt[i:(i + mw - 1)]
    
    if(type=="AR1"){
      nYR <- try(ar.ols(Ytw, aic = FALSE, order.max = 1, demean = TRUE, intercept = FALSE), silent=TRUE)
      if(inherits(nYR,"try-error")) ewsOut[i] <- NA
      else ewsOut[i] <- nYR$ar
    }
    else if(type=="ACF1"){
      nYR <- try(acf(Ytw, lag.max = 1, type = c("correlation"), plot = FALSE), silent=TRUE)
      if(inherits(nYR,"try-error")) ewsOut[i] <- NA
      else ewsOut[i] <- nYR$acf[2]
    }
    else if(type=="densratio"){
      spectfft <- try(spec.ar(Ytw, n.freq = mw, plot = FALSE, order = 1))
      if(inherits(spectfft,"try-error")) ewsOut[i] <- NA
      else ewsOut[i] <- spectfft$spec[2]/spectfft$spec[mw] # low (= 2) vs high (= mw) frequencies
    }
    else if(type=="skew"){
      tmp <- try(abs(moments::skewness(Ytw, na.rm = TRUE)), silent=TRUE)
      if(inherits(tmp,"try-error")) ewsOut[i] <- NA
      else ewsOut[i] <- tmp
    } 
    else if(type=="kurt"){
      tmp <- try(moments::kurtosis(Ytw, na.rm = TRUE), silent=TRUE)
      if(inherits(tmp,"try-error")) ewsOut[i] <- NA
      else ewsOut[i] <- tmp
    } 
    else if(type=="cv"){
      tmp <- try(sd(Ytw, na.rm = TRUE) / mean(Ytw, na.rm = TRUE), silent=TRUE)
      if(inherits(tmp,"try-error")) ewsOut[i] <- NA
      else ewsOut[i] <- tmp
    } 
    else if(type=="sd"){
      tmp <- try(sd(Ytw, na.rm = TRUE), silent=TRUE)
      if(inherits(tmp,"try-error")) ewsOut[i] <- NA
      else ewsOut[i] <- tmp
    } 
    else{
      stop("Invalid type argument")
    }
  }
  
  if(out=="coef"){
    return(ewsOut)
  } 
  else if(out=="trend_kendall"){
    tm<-1:omw
    kdl<-as.numeric(cor.test(tm, ewsOut, alternative="two.sided", method="kendall")$estimate)
    return(kdl)
  }
  else if(out=="trend_spearman"){
    tm<-1:omw
    kdl<-as.numeric(cor.test(tm, ewsOut, alternative="two.sided", method="spearman")$estimate)
    return(kdl)
  }
  else{
    stop("Invalid option in out")
  }

}

#system.time(
#for(i in 1:100) calculate_AR_mw(Yt=rnorm(381),mw_size=30,mw_type="perc",type="skew",out="trend_kendall")
#)



