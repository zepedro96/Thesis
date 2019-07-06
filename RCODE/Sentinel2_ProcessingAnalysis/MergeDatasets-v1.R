
library(raster)
library(readxl)

# Read data ------

Spring16_extrct_v1 <- read.csv("C:/Users/Utilizador/Documents/spring16bands/Spring16_extrct_v1.csv", stringsAsFactors = FALSE)

Passerine <- read.csv("C:/Users/Utilizador/Desktop/TESE/MasterThesisJoseSilva/MasterThesisJoseSilva/DATA/fielddata/FINAL_DATA/PasserineSpeciesDiversity_v1.csv",
                      stringsAsFactors = FALSE)

Flora <- read.csv("C:/Users/Utilizador/Desktop/TESE/MasterThesisJoseSilva/MasterThesisJoseSilva/DATA/fielddata/FINAL_DATA/PlantSpeciesDiversity_Vez_v1.csv",
                  stringsAsFactors = FALSE)

GHC <- read.csv("C:/Users/Utilizador/Desktop/TESE/MasterThesisJoseSilva/MasterThesisJoseSilva/DATA/fielddata/FINAL_DATA/GHC_LifeFormsDiversity_BySSUid.csv")

vezDF <- merge(Passerine, Flora, by="PSU_SSU_ID")
vezDF <- merge(vezDF, Spring16_extrct_v1, by="PSU_SSU_ID") 
vezDF <- vezDF[,-c(40,41,52,53)]

vezDF_GHC <- merge(vezDF, GHC, by="PSU_SSU_ID")

saveRDS(vezDF,"./DATA/vezDF_merge.rds")

vezDF <- readRDS("./DATA/vezDF_merge.rds")

saveRDS(vezDF_GHC,"./DATA/vezDF_GHC_merge.rds")
