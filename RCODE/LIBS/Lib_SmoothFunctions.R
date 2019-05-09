

## Set weights for applying the Whittaker smoother with higher 
## weight for the upper envelope
##
## The weights are set based on annual quantile values
## Slow function to apply over large raster images!!!...

setWeightsByYear<-function(x,indices,quantiles,w,basew=1,check=T,checkWindow=3){
  
  if(length(x) != length(indices))
    stop("Length of x must equal that of indices!")
  
  if(length(quantiles) != length(w))
    stop("Length of quantiles must equal that of weights!")
  
  if(!all(quantiles==sort(quantiles)))
    stop("Quantiles and weights must be ordered in ascendent fashion!")
  
  wts<-rep(basew,length(indices))
  idx<-sort(unique(indices))
  
  for(i in idx){ 
    qts.yr<-quantile(x[indices==i], probs=quantiles)
    
    for(j in 1:length(unique(quantiles))){
      wts[(indices==i) & (x >= qts.yr[j])] <- w[j]   
    }
  }
  
  if(check){
    l<-length(x)
    inc<-floor(checkWindow/2)
    wts_new<-wts
    
    for(i in 1:l){
      win<-c(i-inc,i+inc)
      win<-win[(win>0)&(win<=l)]
      
      if(!any(wts[i]==wts[win])){
        
        med<-mean(wts[win])
        if(med<wts[i]) ## Don't want negative bias to have more weight!..
          wts_new[i]<-med
      }
    }
    return(wts_new)
  }else{
    return(wts)
  }
}



WhittakerAndersonFilter <- function(x, lambda = 1, d = 2){ 
  
  #require(SparseM, warn.conflicts = FALSE)

  #   Smoothing with a finite difference penalty
  #   x:      signal to be smoothed
  #   lambda: smoothing parameter (rough 50..1e4 smooth)
  #   d:      order of differences in penalty (generally 2)
  
  res <- as.numeric(x)

  if(!any(is.na(res))){
    
    m <- length(x)
    E <- as(m, 'matrix.diag.csr')
    class(E) <- 'matrix.csr'
    Dmat <- diff(E, differences = d)
    B <- E + (lambda * t(Dmat) %*% Dmat)
    res <- SparseM::solve(B, res)
    return(res)
    
  }else{
    
      return(NA)
  }
}


hampel.filter <- function(x, k = 7, t0 = 3, ...){
 
  res <- as.vector(x)
  len <- length(x)
  n <- length(which(is.na(res)))
  
  if(n == 0){
    for (i in (k + 1):(len - k)){
      w0 <- res[(i - k):(i + k)]
      x0 <- median(w0, na.rm = TRUE)
      S0 <- 1.4826 * median(abs(w0 - x0), na.rm = TRUE)
      if (abs(res[i] - x0) > (t0 * S0))
        res[i] <- x0
    }
  }
  
  return(res)
}


## Extracts the STL components and returns a vector with all of them
## This function is meant to be used with calc or stackApply from the raster 
## package

stl_extract_components<-function(x, start=c(2000,4), end=c(2017,15), 
                                 frequency=23, s.window=7, na.rm=TRUE,...){
  
  if(any(is.na(x)))
    return(NA)
  
  x.ts<-ts(x,start = start, end=end, frequency=frequency)
  stlObj<-stl(x.ts,s.window = s.window,...)
  
  return(c(as.numeric(stlObj$time.series[,1]), # Seasonal 1 - N
           as.numeric(stlObj$time.series[,2]), # Trend (N+1) : 2*N
           as.numeric(stlObj$time.series[,3])  # Remainder [(2*N)+1] : 3*N
  ))
  
}

# Version with Rcpp Hampel filter
#
#
# smoothTSHampelWhittaker<-function(x, useHampelFilter=TRUE, HF.k=NULL, HF.t0=NULL, HF.winType=NULL, 
#                                   lambda, useWt=TRUE, indices=NULL, quantProbs=NULL, quantWeights=NULL, 
#                                   baseWeight=NULL, checkWts=TRUE, checkWin=3,na.rm=FALSE){
#   
#   if(all(is.na(x))){
#     return(NA)  
#   }else{
#     ## Use Hampel filter/identifier
#     if(useHampelFilter){
#       if(is.null(HF.k) || is.null(HF.t0) || is.null(HF.winType)){
#         stop("Need to set k, t0, and winType for use in Hampel filter/identifier!")
#       }
#     }else{
#       x<-HampelFilterCpp(x, HF.k, HF.t0, HF.winType)
#     } 
#     ## Use weights
#     if(useWt){
#       if(is.null(indices) || is.null(quantProbs) || is.null(quantWeights) || is.null(baseWeight))
#         stop("Parameters indices, quantProbs, quantWeights or baseWeight can not be set to NULL if useWt is TRUE!")
#       ## Rcpp function to weight
#       w<-setWeightsCpp(x,indices, quantProbs, quantWeights, baseWeight, checkWts, checkWin) 
#       return(whit2(x,lambda,w))
#     }else{
#       return(whit2(x,lambda))
#     }
#   }
# }


smoothTS_Hampel_Whittaker<-function(x, useHampelFilter=TRUE, HF.k=NULL, HF.t0=NULL, 
                                  lambda, useWt=TRUE, indices=NULL, quantProbs=NULL, quantWeights=NULL, 
                                  baseWeight=NULL, checkWts=TRUE, checkWin=3,na.rm=FALSE){
  
  if(all(is.na(x))){
    return(NA)  
  }else{
    ## Use Hampel filter/identifier
    if(useHampelFilter){
      if(is.null(HF.k) || is.null(HF.t0)){
        stop("Need to set k and t0 for use in Hampel filter/identifier!")
      }
    }else{
      #x<-HampelFilterCpp(x, HF.k, HF.t0, HF.winType)
      x<-hampel(x, HF.k, HF.t0)$y
    } 
    ## Use weights
    if(useWt){
      if(is.null(indices) || is.null(quantProbs) || is.null(quantWeights) || is.null(baseWeight))
        stop("Parameters indices, quantProbs, quantWeights or baseWeight can not be set to NULL if useWt is TRUE!")
      ## Rcpp function to weight
      w<-setWeightsCpp(x,indices, quantProbs, quantWeights, baseWeight, checkWts, checkWin) 
      return(whit2(x,lambda,w))
    }else{
      return(whit2(x,lambda))
    }
  }
}



