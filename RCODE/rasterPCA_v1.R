library(RStoolbox)
library(raster)
library(sf)
library(ggplot2)

#PCA Raster

shapeArea <- read_sf("D:/Thesis/DATA/Vector/VezBasinCOS12_WGS84_29N.shp")

shapeAreaBuff <- st_buffer(shapeArea, 2500)


fl_10m <- list.files("D:/Thesis/DATA/Raster/SentinelS2A/SENTINEL2A_20160124-114205-331_L2A_T29TNG_D",
                     pattern="_SRE_", full.names = TRUE)[c(3,4,5,9)]

fl_20m <- list.files("D:/Thesis/DATA/Raster/SentinelS2A/SENTINEL2A_20160124-114205-331_L2A_T29TNG_D",
                     pattern="_SRE_", full.names = TRUE)[-c(3,4,5,9)]

rst10m <- stack(fl_10m)
rst20m <- stack(fl_20m)

rst10m <- crop(rst10m, shapeAreaBuff)
rst10m <- mask(rst10m, shapeAreaBuff)

rst20m <- crop(rst20m, shapeAreaBuff)
rst20m <- mask(rst20m, shapeAreaBuff)

rst20m_resamp <- resample(rst20m, rst10m, method="ngb")

rst<- stack(rst10m, rst20m_resamp)

ggRGB(rst, 1,2,3)
set.seed(112)
RstPCA <- rasterPCA(rst,  nSamples = NULL, nComp = 10, spca = FALSE, maskCheck = TRUE,
                    filename = "./OUT/PCA_S2a_Winter2015.tif")

summary(RstPCA$model)

test_PCA <- ggRGB(RstPCA$map,1,2,3, stretch="lin", q=0)



## Make a barplot with explained variance ---------------------------------------- ##

eigs <- RstPCA$model$sdev^2
expVar <- eigs / sum(eigs)
explVarDF <- data.frame(pc = 1:length(expVar),
                        sdev = RstPCA$model$sdev,
                        csum = cumsum(eigs) / sum(eigs),
                        explVar = expVar)

# Proportion of variance explained per principal component
#
g <- ggplot(explVarDF, aes(y=explVar, x=pc)) + 
  geom_bar(stat="identity") + 
  xlab("Principal components") + 
  ylab("Proportion of explained variance")

plot(g)

ggsave("./OUTtoShare/propExplVar_winter15_20161116.png", g)

# Cumulative proportion of explained variance
#
g1 <- ggplot(explVarDF, aes(y=csum, x=pc)) + 
  geom_bar(stat="identity") + 
  xlab("Principal components") + 
  ylab("Cumulative proportion of explained variance")

plot(g1)

ggsave( "./OUTtoShare/cumulPropExplVar_winter15_20161116.png", g1)


