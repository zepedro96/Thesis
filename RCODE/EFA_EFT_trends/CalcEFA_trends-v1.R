

library(raster)
library(fasterize)
library(dplyr)
library(trend)
library(readxl)
library(sf)
library(tidyr)
library(ggplot2)
library(RStoolbox)
library(psych)
library(microbiome)


## ----------------------------------------------------------------------------------------- ##

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

## Functions to calculate diversity and evenness measures -----

divEvFun <- function(x){
  as.data.frame(
    c(evenness(as.numeric(table(x)), index = "all", zeroes = FALSE, detection = 0)[1,,drop=TRUE],
      diversity(as.numeric(table(x)), index = c("shannon","inverse_simpson","gini_simpson"), zeroes = FALSE)[1,,drop=TRUE],
      eft_count = length(unique(x)))
  )
}

evCamargo <- function(x) evenness(as.numeric(table(x)), index = "camargo", zeroes = FALSE, detection = 0)[1,,drop=TRUE]

evSimpson <- function(x) evenness(as.numeric(table(x)), index = "simpson", zeroes = FALSE, detection = 0)[1,,drop=TRUE]

evPielou <- function(x) evenness(as.numeric(table(x)), index = "pielou", zeroes = FALSE, detection = 0)[1,,drop=TRUE]

evEvar<- function(x) evenness(as.numeric(table(x)), index = "evar", zeroes = FALSE, detection = 0)[1,,drop=TRUE]

evBulla <- function(x) evenness(as.numeric(table(x)), index = "bulla", zeroes = FALSE, detection = 0)[1,,drop=TRUE]

divShannon <- function(x) diversity(as.numeric(table(x)), index = "shannon", zeroes = FALSE)[1,,drop=TRUE]

divInvSimpson <- function(x) diversity(as.numeric(table(x)), index = "inverse_simpson", zeroes = FALSE)[1,,drop=TRUE]

divGiniSimpson <- function(x) diversity(as.numeric(table(x)), index = "gini_simpson", zeroes = FALSE)[1,,drop=TRUE]

shannon <- function (x, na.rm=TRUE){
  x <- na.omit(x)
  if(length(x)==0){
    return(0)
  }else{
    tb <- table(x)
    tbrel <- tb / sum(tb)
    z <- tbrel * log(tbrel)
    out <- -1 * sum(z)
    return(out)
  }
}

countDistinct <- function(x){
  x <- na.omit(x)
  if(length(x)==0){
    return(0)
  }else{
    return(length(unique(x))) 
  }
} 

## ----------------------------------------------------------------------------------------- ##

## Functions to calculate the Theil-Sen trend slope and the p-value

senSlope <- function(x,...) sens.slope(x,...)$estimates
senSlopePval <- function(x,...) sens.slope(x,...)$p.value


## ----------------------------------------------------------------------------------------- ##

## Load UTM 1km2 grid -----

#grid1km   <- read_sf("D:/Dropbox/MasterThesisJoséSilva/DATA/Grids/SampGridVez_1km_Vez_All_v1.shp")
grid1km   <- read_sf("D:/MyDocs/Dropbox/MasterThesisJoséSilva/DATA/Grids/SampGridVez_1km_Vez_All_v1.shp")



## ----------------------------------------------------------------------------------------- ##
## EFA's -----
## ----------------------------------------------------------------------------------------- ##

# List geotiff files for EFT's by year
#
fl <- list.files("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun", pattern=".tif$", full.names = TRUE)
bn <- basename(fl)

# Names for each EFA variable
varNames <- unlist(lapply(strsplit(bn,"_"), FUN = function(x) gsub(".tif","",x[length(x)])))

# Short names for each EFA variable (must length(varNames) == length(smallNames) and 
# use the same order exactly)
smallNames <- c("am","av","cv","dq","mx","md","mi","mm","sx","sp","sd","wt","ei")
names(smallNames) <- varNames

outFolder <- "./OUTtoSHARE/EFA_trends"


# Define ranges of variables for making the trend analysis in three nested periods
# One entry in the list for each measure used for spatial aggregation/diversity calculation
#
#
colRanges <- list(
  # Spatial standard-deviation by grid unit
  "spSTD" = list(
  "all"             = c(3:36),
  "last_ten_yrs"    = c(27:36),
  "last_twenty_yrs" = c(17:36)
  ), 
  # Spatial average by grid unit
  "spAVG" = list(
    "all"             = c(39:72),
    "last_ten_yrs"    = c(63:72),
    "last_twenty_yrs" = c(53:72)
  )
)

pb <- txtProgressBar(1,length(fl),style=3)


# Iterate by annual EFA
#
#
for(i in 1:length(fl)){
  
  vname <- varNames[i]
  sname <- smallNames[vname]
  
  dir.create(paste(outFolder,vname,sep="/"))
  
  EFA_rst <- stack(fl[i])
  names(EFA_rst) <- paste(sname,substr(as.character(1984:2019),3,4),sep="_")
  
  birdGridRst1km <- fasterize::fasterize(grid1km, EFA_rst[[1]], field = "ID_PSU")
  
  rstStack <- stack(birdGridRst1km, EFA_rst)
  
  rstDF <- values(rstStack) %>% na.omit
  colnames(rstDF)[1] <- "ID_PSU"
  
  # Aggregate data from pixel-level to 1km2 units
  rstDF_agg <- rstDF %>%
    as.data.frame %>%
    group_by(ID_PSU) %>%
    summarize_all(.funs=list("sd" = sd, 
                             "av" = mean))
  
  # Iterate by spatial aggregation metric (e.g., sd, mean, etc...)
  #
  #
  for(j in 1:length(colRanges)){
    
    spVarName <- names(colRanges)[j]
    
    dir.create(paste(outFolder,vname,spVarName,sep="/"))
    
    ## EFA's Estimate trends over different periods -------------------------------------------- ##
    
    r_all <- colRanges[[j]][["all"]]
    r_10y <- colRanges[[j]][["last_ten_yrs"]]
    r_20y <- colRanges[[j]][["last_twenty_yrs"]]
    
    # print(colnames(rstDF_agg)[r_all])
    # print(colnames(rstDF_agg)[r_10y])
    # print(colnames(rstDF_agg)[r_20y])
    
    # Trend for the whole analysis period
    trend_all <- rstDF_agg %>% 
      select(r_all) %>% 
      data.frame(sen_all = apply(.,1,senSlope),pv_all = apply(.,1,senSlopePval)) %>% 
      select(sen_all, pv_all)
    
    # Last ten years trend
    trend_10yr <- rstDF_agg %>% 
      select(r_10y) %>% 
      data.frame(sen_10yr = apply(.,1,senSlope),pv_10yr = apply(.,1,senSlopePval)) %>% 
      select(sen_10yr, pv_10yr)
    
    # Last twenty years trend
    trend_20yr <- rstDF_agg %>% 
      select(r_20y) %>% 
      data.frame(sen_20yr = apply(.,1,senSlope),pv_20yr = apply(.,1,senSlopePval)) %>% 
      select(sen_20yr, pv_20yr)
    
    # Combine all the results
    trendsDF <- data.frame(ID_PSU=rstDF_agg$ID_PSU, trend_10yr, trend_20yr, trend_all)
    
    # Append all the results to the vector grid
    grid1km_append <- grid1km %>% 
      left_join(rstDF_agg, by = "ID_PSU")%>% 
      left_join(trendsDF, by = "ID_PSU")
    
    outFn1 <- paste(outFolder,"/",vname,"/",spVarName,"/EVI_",vname,"_",spVarName,"_byYear_grid1km.shp",sep="")
    outFn2 <- paste(outFolder,"/",vname,"/",spVarName,"/EVI_",vname,"_",spVarName,"_byYear_grid1km_centrds.shp",sep="")
    
    write_sf(grid1km_append,outFn1)
    write_sf(grid1km_append %>% st_centroid,outFn2)

  }
  
  setTxtProgressBar(pb, i)
}



## ----------------------------------------------------------------------------------------- ##
## EFT's -----
## ----------------------------------------------------------------------------------------- ##


# List geotiff files for EFT's by year
#
fl <- list.files("./DATA/RASTER/Landsat/EVI-32day/metrics/jul-jun/EFT", pattern=".tif$", full.names = TRUE)
EFA_rst <- stack(fl)

outFolder <- "./OUTtoSHARE/EFT_trends"


EFA_name <- "et" ## USE TWO LETTERS ONLY (shapefiles will not work with more characters... :-0
names(EFA_rst) <- paste(EFA_name,substr(as.character(1984:2019),3,4),sep="_")

birdGridRst1km <- fasterize::fasterize(grid1km, EFA_rst[[1]], field = "ID_PSU")

rstStack <- stack(birdGridRst1km, EFA_rst)

rstDF <- values(rstStack) %>% na.omit
colnames(rstDF)[1] <- "ID_PSU"

# Aggregate data from pixel-level to 1km2 units
# Calculate evenness and diversity metrics by grid unit
rstDF_agg <- rstDF %>% 
  as.data.frame %>% 
  group_by(ID_PSU) %>% 
  summarize_all(.funs=list(# Evenness indices
                           "ec" = evCamargo,
                           "es" = evSimpson,
                           "ep" = evPielou,
                           "ee" = evEvar,
                           "eb" = evBulla,
                           # Diversity indices
                           "ds" = divShannon,
                           "di" = divInvSimpson,
                           "dg" = divGiniSimpson,
                           "rc" = countDistinct))


## EFT's Estimate trends over different periods -------------------------------------------- ##



# Define names for each variable/evenness/diversity measure
varNames <- c("evCamargo","evSimpson","evPielou","evEvar", "evBulla",
            "divShannon","divInvSimpson","divGiniSimpson","EFTcount")


# Define column-wise ranges for each evenneness and diversity index
chunks <- chunk2(2:ncol(rstDF_agg), length(varNames)) # 9 variables in total

i <- 0
for(chunk in chunks){
  
  i <- i + 1
  v <- varNames[i]
  
  dir.create(paste(outFolder,v,sep="/"))
  
  cat("Doing chunk",i,"out of",length(chunks),"\n\n")
  
  # Define ranges for each analysis period: 10-yr, 20-yr and all 
  # (NOTE: avoid start and end incomplete years)
  r_10y <- chunk[26:35]
  r_20y <- chunk[16:35]
  r_all <- chunk[2:35]
  
  print(colnames(rstDF_agg)[r_10y])
  cat("\n\n")
  
  print(colnames(rstDF_agg)[r_20y])
  cat("\n\n")
  
  print(colnames(rstDF_agg)[r_all])
  cat("\n\n")
  
  # Calculate trends for the three nested periods
  # Last ten years trend
  trend_10yr <- rstDF_agg %>% 
    select(r_10y) %>% 
    data.frame(sen_10yr = apply(.,1,senSlope),pv_10yr = apply(.,1,senSlopePval)) %>% 
    select(sen_10yr, pv_10yr)
  
  # Last twenty years trend
  trend_20yr <- rstDF_agg %>% 
    select(r_20y) %>% 
    data.frame(sen_20yr = apply(.,1,senSlope),pv_20yr = apply(.,1,senSlopePval)) %>% 
    select(sen_20yr, pv_20yr)
  
  # All years trend
  trend_all <- rstDF_agg %>% 
    select(r_all) %>% 
    data.frame(sen_all = apply(.,1,senSlope),pv_all = apply(.,1,senSlopePval)) %>% 
    select(sen_all, pv_all)
  
  # Append results for all trends estimated
  trendsDF <- data.frame(ID_PSU=rstDF_agg$ID_PSU, trend_10yr, trend_20yr, trend_all)
  
  # Append trend data to the vector grid
  grid1km_append <- grid1km %>% 
    left_join(rstDF_agg[, c(1,chunk)], by = "ID_PSU")%>% 
    left_join(trendsDF, by = "ID_PSU")
  
  # Write data files
  outFn1 <- paste(outFolder,"/",v,"/EFT_",varNames[i],"_grid1km_1984_2019_v1.shp",sep="")
  outFn2 <- paste(outFolder,"/",v,"/EFT_",varNames[i],"_grid1kmCentrs_1984_2019_v1.shp",sep="")
  
  write_sf(grid1km_append, outFn1)
  write_sf(grid1km_append %>% st_centroid, outFn2)

}


