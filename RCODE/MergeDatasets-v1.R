
library(raster)
library(readxl)

# Extract data from raster image ------

extRst_SD <- extract(rst, birdsgrid, fun=sd)
colnames(extRst_SD) <- paste(colnames(extRst_SD),"_sd",sep="")
colnames(extRst_AVG) <- paste(colnames(extRst_AVG),"_avg",sep="")
extRst <- cbind(birdsgrid@data[,c("ID_PSU","ID_SSU","PSU_SSU_ID")], extRst_AVG, extRst_SD)
write.csv(extRst, "Spring16_extrctBands.csv",row.names = FALSE)

# Read excel data ------

read_excel("C:/Users/Utilizador/Desktop/TESE/Extracts.xlsx")
extNDVI <- read_xlsx("C:/Users/Utilizador/Desktop/TESE/Extracts_NDVI.xlsx", sheet = "mean_&_stdev_NDVI", col_names = TRUE)
TotExt <- merge.data.frame(extRst, extNDVI)
write.csv(TotExt, "Spring16_extrct.csv",row.names = FALSE)

Spring16_extrct_v1 <- read.csv("C:/Users/Utilizador/Documents/spring16bands/Spring16_extrct_v1.csv", stringsAsFactors = FALSE)

Passerine <- read.csv("C:/Users/Utilizador/Desktop/TESE/MasterThesisJoseSilva/MasterThesisJoseSilva/DATA/fielddata/FINAL_DATA/PasserineSpeciesDiversity_v1.csv",
                      , stringsAsFactors = FALSE)

Flora <- read.csv("C:/Users/Utilizador/Desktop/TESE/MasterThesisJoseSilva/MasterThesisJoseSilva/DATA/fielddata/FINAL_DATA/PlantSpeciesDiversity_Vez_v1.csv",
                  , stringsAsFactors = FALSE)


vezDF <- merge(Passerine, Flora, by="PSU_SSU_ID")
vezDF <- merge(vezDF, Spring16_extrct_v1, by="PSU_SSU_ID") 
vezDF <- vezDF[,-c(40,41,52,53)]
