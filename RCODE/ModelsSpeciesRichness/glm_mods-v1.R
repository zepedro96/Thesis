

library(DescTools)
library(car)

vezDF_vars1km <- read.csv("./OUT/vezDF_vars1km_v2.csv")

varsToUse <- c("EFAavg_std","EFAsemax_std","evar","eft_count")

#combns<- combn(1:4,2)

respVars <- colnames(vezDF_vars1km)[c(3, 7:13, 15:32)] # Select response variables
#respVar <- "Sp_rich_sum"

outMods <- list()

results <- as.data.frame(matrix(NA,ncol=4, nrow = length(respVars)))
colnames(results) <- c("respVar","CoxSnell","Nagelkerke","Effron")


vifRes <- as.data.frame(matrix(NA,ncol=length(varsToUse), nrow = length(respVars)))
colnames(vifRes) <- varsToUse


for(i in 1:length(respVars)){
  
  # v1<-combns[1,i]
  # v2<-combns[2,i]
  # var1<- varsToUse[v1]
  # var2<- varsToUse[v2]
  
  respVar <- respVars[i]
  #form <- paste(respVar,paste(var1,var2,sep="+"),sep="~")
  form <- paste(respVar,paste(varsToUse,collapse="+"),sep="~")
  
  compForm <- as.formula(form)
  
  summary(form)
  
  mod <- glm(formula = compForm, data = vezDF_vars1km, family = poisson())
  
  vifRes[i,] <- sqrt(vif(mod))
  #print(mod)
  
  outMods[[i]] <- mod
  
  results[i,1] <- respVar
  results[i, 2:4] <- PseudoR2(mod,"all")[c(3,4,7)]
  
  cat("Finished model:",i,"\n\n")
  
}

write.csv(vifRes, "./OUT/vifRes_v1.csv")
write.csv(results, "./OUT/results_glm_v2.csv")
