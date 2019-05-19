
library(raster)
library(dplyr)

LT5 <- stack("./DATA/RASTER/Landsat/EVI-32day/LT5_EVI_Vez_Basin_1984-2012.tif")
LT7 <- stack("./DATA/RASTER/Landsat/EVI-32day/LT7_EVI_Vez_Basin_1999-2019.tif")
LT8 <- stack("./DATA/RASTER/Landsat/EVI-32day/LT8_EVI_Vez_Basin_2013-2019.tif")

LTdates <- readxl::read_excel("./DATA/RASTER/Landsat/EVI-32day/EVI-32-day_ImageList-v1.xlsx")

LT5_dates <- LTdates %>% filter(Satellite=="LT05")
LT7_dates <- LTdates %>% filter(Satellite=="LE07")
LT8_dates <- LTdates %>% filter(Satellite=="LC08")

LTdates %>% 
  group_by(Satellite, Year) %>% 
  summarize(nimgs=n())

LT5_dates %>% 
  group_by(Satellite, Year) %>% 
  summarize(nimgs=n()) %>% 
  as.data.frame()

DT_all <- rbind(
  #LT5_dates[[1:283]],
  LT5_dates[1:220,],
  LT7_dates[39:49,],
  LT5_dates[232:283,],
  LT7_dates[101:102,],
  LT5_dates[284:331,],
  LT7_dates[151:154,],
  LT5_dates[332,],
  LT7_dates[156:166,],
  LT8_dates[,]
  )

DT_all %>% group_by(Year) %>% summarize(nimgs=n()) %>% as.data.frame()

DT_all <- DT_all %>% mutate(fileName = paste("EVI_T1_32day_",DateCode,"_",Satellite,"-v1.tif",sep=""))

#View(DT_all)
write.csv(DT_all, "./DATA/RASTER/Landsat/EVI-32day/FullImageList_1984-2019-LT-5-7-8_v2.csv",row.names = FALSE)

LT_all <- stack(
  #LT5[[1:283]],
  LT5[[1:220]],
  LT7[[39:49]],
  LT5[[232:283]],
  LT7[[101:102]],
  LT5[[284:331]],
  LT7[[151:154]],
  LT5[[332]],
  LT7[[156:166]],
  LT8
)

print(nrow(DT_all))
print(nlayers(LT_all))



## Check if the files in the folder are in the list??!

fl <- list.files("./DATA/RASTER/Landsat/EVI-32day/full",pattern=".tif$")

wdir <- "./DATA/RASTER/Landsat/EVI-32day/full/"



for(i in 1:length(fl)){
  
  fn <- fl[i]
  
  if(!(fn %in% DT_all$fileName)){
    cat("\n\n->",fn,"not in the base file list!! Removing it....")
    try(file.remove(paste(wdir,fn,sep="")))
    cat(" done!\n")
  }
}

#writeRaster(LT_all,"./DATA/RASTER/Landsat/EVI-32day/EVI_32day_VezBasin_1984_2019_LT578-v1.tif")

N <- nlayers(LT_all)
pb <- txtProgressBar(1,N,style = 3)

for(i in 1:N){
  
  dateCode <- DT_all$DateCode[i]
  satCode <- DT_all$Satellite[i]
  
  outFileName <- paste(wdir, "EVI_T1_32day_",dateCode,"_",satCode,"-v1.tif",sep="")
  
  if(!file.exists(outFileName)){
    cat("\n-> Writing file:",outFileName,"\n\n")
    writeRaster(LT_all[[i]],outFileName)
  }else{
    cat("\n-> Skipping file!!\n\n")
  }
  setTxtProgressBar(pb, i)
}

