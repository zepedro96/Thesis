
library(dplyr)

birdDF <- read.csv("C:/Users/JG/Desktop/birdData_GBIF/dados_gbif_vez_v2.txt")

countDistinct <- function(x) length(na.omit(unique(x))) 


vals <- data.frame(minVs=1:250, maxVs=c((1:250)+5, (1:250)+8, (1:250)+10, 
                                         (1:250)+12, (1:250)+15, (1:250)+17,
                                         (1:250)+20, (1:250)+30, (1:250)+50))

res <- vector(mode="numeric",length=nrow(vals))
nr  <- vector(mode="numeric",length=nrow(vals))


for(i in 1:nrow(vals)){
  
  DF <- try(birdDF %>% 
    group_by(ID_PSU) %>% 
    summarize(SpRich = countDistinct(species), npts=n()) %>% 
    filter(npts >= vals[i,1], npts <= vals[i,2])) 
  
    if(nrow(DF) >= 10){
    
      corVal<- DF %>% 
        cor %>% round(2) %>% `[`(3,2)
      corVal<-ifelse(inherits(corVal,"try-error"),NA,corVal)
      res[i]<-corVal
      nr[i]<-nrow(DF)
  } 
}

data.frame(corVal=res, nr=nr, vals) %>% 
  arrange(desc(nr)) %>% 
  filter(corVal>0, corVal<0.5)
  
DF <- try(birdDF %>% 
            group_by(ID_PSU) %>% 
            summarize(SpRich = countDistinct(species), npts=n()) %>% 
            filter(npts >= 20, npts <= 37)) %>% 
            arrange(desc(SpRich))

write.csv(DF, "C:/Users/JG/Desktop/birdData_GBIF/gbif_data_sprich_v2.csv", row.names = FALSE)

