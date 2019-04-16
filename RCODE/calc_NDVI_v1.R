#median, MaxVal & mad

fileList <- list.files(path="DATA/Raster/S2A_crop",  full.names=TRUE, recursive=FALSE, pattern = "_NDVI")

stackNDVI <- stack(fileList)

median <- calc(stackNDVI, fun = median)

MaxVal <- calc(stackNDVI, fun = maxValue

mad <- mad(stackNDVI, center = median, constant = 1.4826)
