

library(raster)
library(sf)

shapeArea <- read_sf("shape_path.shp")

shapeAreaBuff <- st_buffer(shapeArea, 5000)

rst <- stack("file")
rst <- crop(rst, shapeAreaBuff)
rst <- mask(rst, shapeAreaBuff)

ndvi <- (rst[[7]] - rst[[3]]) / (rst[[7]] + rst[[3]])

writeRaster(ndvi, "filename.tif")