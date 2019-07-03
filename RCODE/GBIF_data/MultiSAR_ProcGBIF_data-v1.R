

library(dplyr)
library(sf)
library(sp)
library(mmSAR)

#rm(list=ls())

data(power)
data(expo)
data(negexpo)
data(monod)
data(ratio)
data(logist)
data(lomolino)
data(weibull)


pts <- read_sf("C:/Users/JG/Desktop/GBIF_DATA/GBIF_Passeriformes_NW_PTES_29N_v1.shp")

# Target areas
tareas <- seq(1, 25, by=1)

# Radius of circunference in meters
rs <- sqrt((tareas / pi))*1000

maxDist <- units::set_units(max(rs),m)

rids <- 1:nrow(pts)

coords <- st_coordinates(pts)

xmin <- min(coords[,1])
xmax <- max(coords[,1])

ymin <- min(coords[,2])
ymax <- max(coords[,2])



nr <- 0
SARobjs <- list()

n_find <- 20
n_max <- 2000

it <- 0

mods <- c("power","expo","negexpo","monod","logist",
          "ratio","lomolino","weibull")

repeat{
  
  it <- it + 1
  xs <- runif(1, xmin, xmax)
  ys <- runif(1, ymin, ymax)
  
  pt <- st_sfc(st_point(c(xs,ys)))
  st_crs(pt) <- st_crs(pts)
  
  # Stage 1 selection: 
  #   1) Randomly place point, then draw a 5 km radius
  #   2) In that 5km radius, pick a random point  (if any)
  #
  
  selIDs_st1 <- st_intersects(pts, st_buffer(pt, 5000), sparse = FALSE)
  
  nspMat <- matrix(NA,nrow = length(rs), ncol=3)
  
  if(sum(selIDs_st1[,1]) > 0){
    
    selID_st2 <- sample(rids[selIDs_st1[,1]],1)
    selPoint_st2 <- pts[selID_st2,]
    
    # Check if the selected coordinates overlap their buffer areas
    #
    if(exists("centralCoords")){
      
      centralPoints_sf <- st_as_sfc(SpatialPoints(centralCoords))
      st_crs(centralPoints_sf) <- st_crs(pts)
      
      dists <- st_distance(selPoint_st2, centralPoints_sf)[1,]
      
      if(min(dists) < 2 * maxDist){
        cat("Skipping to next iteration - minimum distance condition not verified!!\n\n")
        next
      }
    }

    pb <- txtProgressBar(1, length(rs), style = 3)
    
    for(i in 1:length(rs)){
      
      buffRadius <- rs[i]
      
      selIDsInRange <- st_intersects(pts, st_buffer(selPoint_st2, buffRadius), sparse = FALSE)
      ptsInRange <- pts[selIDsInRange[,1], ]
      
      species <- ptsInRange %>% 
        select(species) %>% 
        st_drop_geometry %>% 
        pull %>% 
        unique %>% 
        na.omit
      
      nspMat[i,] <- c(rs[i],tareas[i],length(species))
      
      setTxtProgressBar(pb, i)
    }
    
    print(nspMat)
    print(sd(nspMat[,3]))
    print(min(nspMat[,3]))
    
    cofv <- (sd(nspMat[,3])/mean(nspMat[,3])) * 100
    print(cofv)
    
    if(cofv >= 20 & min(nspMat[,3]) < 30){
      
      
      sarData <- list(name = paste("PLOT",nr,sep="_"), data = data.frame(a = nspMat[,2], s = nspMat[,3]))
      
      mtSAR <- try(multiSAR(mods, sarData, nBoot=99, crit="Info", norTest="lillie", verb=FALSE))
      
      if(!inherits(mtSAR,"try-error")){
        
        nr <- nr+1

        plot(nspMat[,2],nspMat[,3])
        lines(nspMat[,2],mtSAR$averaged)
        
        SARobjs[[nr]] <- mtSAR
        tmpSpRichSAR <- data.frame(nr = nr, a = nspMat[,2], s = nspMat[,3], s_avg = mtSAR$averaged)
        tmpCoords <- st_coordinates(selPoint_st2)
        
        if(nr==1){
          centralCoords <- tmpCoords
          spRichSAR <- tmpSpRichSAR
        }else{
          centralCoords <- rbind(centralCoords,tmpCoords)
          rownames(centralCoords) <- 1:nrow(centralCoords)
          
          spRichSAR <- rbind(spRichSAR, tmpSpRichSAR)
        }
        
        cat("\n\n-> Found point nr [",nr,"]\n\n")
        
      }
    }
  }
  
  cat("\n\n### Finished round nr:",it,"[n_found =",nr,"]... ###\n\n")
  
  # Check stopping conditions
  if(nr==n_find | it == n_max)
    break
}

write.csv(centralCoords, "./DATA/SAR_CentralPointCoords-v1.csv", row.names = FALSE)
write.csv(spRichSAR, "./DATA/SAR_spRichObsVsEstimated-v1.csv", row.names = FALSE)

save.image(file="./DATA/SAR_analysis-v1.RData")

spDF <- st_as_sf(SpatialPointsDataFrame(centralCoords, data = as.data.frame(centralCoords)))
st_crs(spDF) <- st_crs(pts)

write_sf(spDF,"./DATA/SAR_CentralPointCoords_v1.shp")

