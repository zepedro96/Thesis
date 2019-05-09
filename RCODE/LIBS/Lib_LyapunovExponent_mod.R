

## João Gonçalves
## Porto, 2016, 2017
##
## Ecosystem stability analyses of Portugal protected areas
##
## Changes the original functions from tseriesChaos package to avoid 
## messages in raster calc operations
##

checkEmbParms <- function(series, m, d, t=0, s=1, ref=NULL) {
  n <- length(series)-(m-1)*d
  if(n <= 0) stop("Not enough points to handle these parameters")
  if(!is.null(ref)) if (ref>n) stop("Not enough points to handle these parameters")
  if(t<0) stop("Theiler window t must be non-negative")
  if(s<=0) stop("Number of steps must be positive")
  
}

find_knearests <- function(series, m, d, t, eps, ref, k, s){
  x <- ( series - min(series) ) / (diff(range(series)))
  eps <- eps/diff(range(series))
  res <- numeric(ref*k)
  res <- .C("find_knearests", series=as.double(x), m=as.integer(m), d=as.integer(d), t=as.integer(t), 
            length=as.integer(length(x)), eps=as.double(eps),ref=as.integer(ref), k=as.integer(k), 
            s=as.integer(s), res=as.integer(res), PACKAGE="tseriesChaos")[["res"]]
  res[res==-1] <- NA;
  matrix(res, ref, k)
}

follow_points <- function(series, m, d, ref, k, s, nearest) {
  res <- numeric(s)
  nearest[is.na(nearest)] <- -1
  .C("follow_points",series=as.double(series), m=as.integer(m), d=as.integer(d), length=as.integer(length(series)), 
     nref=as.integer(length(ref)), nrow=as.integer(nrow(nearest)), k=as.integer(k), s=as.integer(s), 
     nearest=as.integer(nearest), ref=as.integer(ref), res=as.double(res),PACKAGE="tseriesChaos")[["res"]]
}

## Removes messages from the function to increase processing speed
##
LyapunovExp<-function (series, m, d, t, k = 1, ref, s, eps){
  
  checkEmbParms(series, m, d, t, s, ref)
  series <- as.ts(series)
  n <- length(series) - (m - 1) * d - s
  
  if (ref < 0) 
    ref <- n
  trash <- numeric()
  ref <- 1:ref

  nearest <- find_knearests(series, m = m, d = d, t = t, ref = length(ref), s = s, eps = eps, k = k)
  trash <- apply(nearest, 1, function(x) any(is.na(x)))
  ref <- ref[!trash]
  if (length(ref) == 0) 
    stop("Not enough neighbours found")

  res <- follow_points(series, m = m, d = d, s = s, ref = ref, nearest = nearest, k = k)
  ts(res, frequency = frequency(series), start = 0)
}

