

library(DescTools)
library(car)
library(MASS)
library(glmnet)
library(DCluster)

#vezDF_vars1km <- read.csv("./OUT/Landsat/vezDF_vars1km_v2.csv")
vezDF_vars1km <- read.csv("./OUT/vezDF_vars1km_v2.csv")

varsToUse <- c(
               "EFAmax_std","EFAamp_std","EFAavg_std","EFAmin_std" ,"EFAsemax_std",
               "EFAstd_std",
               "evar","shannon","eft_count",
               "EFAmin_avg" ,"EFAmax_avg" ,"EFAavg_avg","EFAsemax_avg",
               "EFAstd_avg","EFAamp_avg"
               )

#varsToUse <- colnames(vezDF_vars1km)[32:66]

respVars <- colnames(vezDF_vars1km)[c(3, 7:13, 15:31)] # Select response variables
#respVar <- "Sp_rich_sum"

outMods <- list()

results <- as.data.frame(matrix(NA,ncol=7, nrow = length(respVars)))
colnames(results) <- c("respVar","CoxSnell","Nagelkerke","Effron",
                       "OdT.LRT","OdT.DeanB","OdT.DeanB2")

vifRes <- as.data.frame(matrix(NA,ncol=length(varsToUse), 
                               nrow = length(respVars)))
colnames(vifRes) <- varsToUse
rownames(vifRes) <- respVars

nsplits <- 5
nvarsMax <- 4


for(i in 1:length(respVars)){
  
  respVar <- respVars[i]
  
  glmLassoReg <- glmnet(x= scale(as.matrix(vezDF_vars1km[,varsToUse]), center = T, scale = T),
                        y=vezDF_vars1km[,respVar], standardize = FALSE,
                        family = "poisson", alpha = 0.5, dfmax = nvarsMax)
  
  tmpLasso <- apply(as.matrix(glmLassoReg$beta),1,max)
  
  # cvGLMnet <- cv.glmnet(x=as.matrix(vezDF_vars1km[,varsToUse]),
  #           y=vezDF_vars1km[,respVar], 
  #           family = "poisson", alpha = 0.5, nfolds=5)
  # predict(cvGLMnet, newx = as.matrix(vezDF_vars1km[,varsToUse]), type="response")
  
  
  splits <- split(1:nrow(vezDF_vars1km),f = 1:nsplits)
  
  cvRes <- as.data.frame(matrix(NA,length(splits), ncol=4))
  
  for(z in 1:10){
    for(k in 1:length(splits)){
      
      inds <- splits[[k]]
      
      glmLassoRegCV <- glmnet(x = scale(as.matrix(vezDF_vars1km[-inds,varsToUse]), center = TRUE, scale = TRUE),
                            y = vezDF_vars1km[-inds,respVar], dfmax = nvarsMax,
                            family = "poisson", alpha = 0.5, standardize = FALSE)
      
      pred <- predict(glmLassoRegCV, newx = as.matrix(vezDF_vars1km[inds,varsToUse]),type="response")
      
      obs <- vezDF_vars1km[inds,respVar]
      
      R2.fold <- cor(obs,pred[,ncol(pred)])^2
      Scor.fold <- cor(obs,pred[,ncol(pred)], method="spearman")
      
      cvRes[k,1] <- respVar
      cvRes[k,2:4] <- c(k, R2.fold, Scor.fold)
    }
    
    if(z==1){
      cvRes1 <- cvRes
    }else{
      cvRes1 <- rbind(cvRes1, cvRes)
    }
  }
  
  if(i==1){
    lassoCoeffs <- tmpLasso
    cvResFinal <- cvRes1
  }else{
    cvResFinal <- rbind(cvResFinal,cvRes1)
    lassoCoeffs <- rbind(lassoCoeffs,tmpLasso)
  }
  
  nonNullCoeffs <- sort((tmpLasso[tmpLasso > 0]), decreasing = TRUE)
  
  if(length(nonNullCoeffs)!=0){
    if(length(nonNullCoeffs) < nvarsMax){
      lassoVars <- names(nonNullCoeffs[1:length(nonNullCoeffs)])
    }else{
      lassoVars <- names(nonNullCoeffs[1:nvarsMax])
    }
  }else{
    next
  }
  
  form <- paste(respVar,paste(lassoVars, collapse="+"),sep="~")
  compForm <- as.formula(form)
  
  modDF <- data.frame(vezDF_vars1km[,respVar,drop=FALSE], 
                      scale(vezDF_vars1km[,lassoVars,drop=FALSE], center=TRUE,scale=TRUE))
  
  mod <- glm(formula = compForm, data = modDF, family = poisson())

  if(length(lassoVars) > 2){
    vifRes[i, lassoVars] <- vif(mod)
  }else{
    vifRes[i, lassoVars] <- NA
  }
  
  coeffs <- summary(mod)$coefficients
  tmp <- data.frame(respVar=respVar, coeffs, sig=as.integer(coeffs[,4] < 0.05))
  
  if(i==1){
    modSummary <- tmp
  }else{
    modSummary <- rbind(modSummary,tmp)
  }
  
  mod.nb <- try(glm.nb(compForm, data = vezDF_vars1km))
  if(!inherits(mod.nb,"try-error")){
    test0 <- test.nb.pois(mod.nb, mod)
    test0 <- test0$p.value
  }else{
    test0 <- NA
  }
  
  test1 <- DeanB(mod)
  test2 <- DeanB2(mod)
  
  outMods[[i]] <- mod
  
  results[i,1] <- respVar
  results[i, 2:4] <- PseudoR2(mod,"all")[c(3,4,7)]
  results[i, 5:7] <- c(test0, test1$p.value, test2$p.value)
  
  print(respVar)
  print(nonNullCoeffs)
  cat("Finished model:",i,"\n\n")
  
}

colnames(lassoCoeffs) <- varsToUse
rownames(lassoCoeffs) <- respVars
#lassoCoeffs <- round(lassoCoeffs,5)
#lassoCoeffs <- rbind(lassoCoeffs,Count=apply(lassoCoeffs,2,FUN = function(x) sum(x!=0)))


modSelFreq <- sort(apply(lassoCoeffs,2,FUN = function(x) sum(x!=0)), decreasing = TRUE)
modSelFreq <- data.frame(varName = names(modSelFreq), freq = modSelFreq, relfreq=(modSelFreq/nrow(lassoCoeffs))*100)
print(modSelFreq)


write.csv(modSelFreq, "./OUTtoShare/varSelectionFrequency_elasticnet_v5.csv",row.names = FALSE)

write.csv(modSummary, "./OUTtoShare/ModelsSummaries-top3vars_Poisson-GLM_v5.csv")
write.csv(lassoCoeffs, "./OUTtoShare/ModelSelection_CoeffsElasticnet_v5.csv")

write.csv(vifRes, "./OUTtoShare/VIF_Results-top3vars_glm_elasticnet_v5.csv")
write.csv(results, "./OUTtoShare/Results_GLM_elasticnet_v5.csv", row.names = FALSE)


