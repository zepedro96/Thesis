

library(DescTools)
library(car)
library(MASS)
library(glmnet)
library(DCluster)

vezDF_vars1km <- read.csv("./OUT/Landsat/vezDF_vars1km_v2.csv")

varsToUse <- c("EFAmin_std" ,"EFAmax_std" ,"EFAavg_std","EFAsemax_std",
               "EFAstd_std","EFAamp_std",
               
               "EFAmin_avg" ,"EFAmax_avg" ,"EFAavg_avg","EFAsemax_avg",
               "EFAstd_avg","EFAamp_avg",
               
               "evar","shannon","eft_count")

#varsToUse <- colnames(vezDF_vars1km)[32:66]

respVars <- colnames(vezDF_vars1km)[c(3, 7:13, 15:32)] # Select response variables
#respVar <- "Sp_rich_sum"

outMods <- list()

results <- as.data.frame(matrix(NA,ncol=7, nrow = length(respVars)))
colnames(results) <- c("respVar","CoxSnell","Nagelkerke","Effron",
                       "OdT.LRT","OdT.DeanB","OdT.DeanB2")



vifRes <- as.data.frame(matrix(NA,ncol=length(varsToUse), 
                               nrow = length(respVars)))
colnames(vifRes) <- varsToUse
rownames(vifRes) <- respVars



for(i in 1:length(respVars)){
  
  respVar <- respVars[i]
  # form <- paste(respVar,paste(varsToUse,collapse="+"),sep="~")
  # compForm <- as.formula(form)
  
  glmLassoReg <- glmnet(x=as.matrix(vezDF_vars1km[,varsToUse]),
                        y=vezDF_vars1km[,respVar], 
                        family = "poisson", dfmax=3)
  
  tmpLasso <- apply(as.matrix(glmLassoReg$beta),1,max)
  if(i==1){
    lassoCoeffs <- tmpLasso
  }else{
    lassoCoeffs <- rbind(lassoCoeffs,tmpLasso)
  }
  
  nonNullCoeffs <- (tmpLasso[tmpLasso > 0])
  
  if(length(nonNullCoeffs)!=0){
    if(length(nonNullCoeffs)<3){
      lassoVars <- names(nonNullCoeffs[1:length(nonNullCoeffs)])
    }else{
      lassoVars <- names(nonNullCoeffs[1:3])
    }
  }else{
    next
  }
  
  form <- paste(respVar,paste(lassoVars, collapse="+"),sep="~")
  compForm <- as.formula(form)
  
  mod <- glm(formula = compForm, data = vezDF_vars1km, family = poisson())
  #mod.step <- stepAIC(mod, direction = "both")

  if(length(lassoVars)>2){
    vifRes[i, lassoVars] <- vif(mod)
  }else{
    vifRes[i, lassoVars] <- NA
  }
  
  
  coeffs <- summary(mod)$coefficients
  tmp <- data.frame(respVar=respVar, coeffs, sig=as.integer(coeffs[,4] < 0.05))
  
  if(i==1){
    modSummary<-tmp
  }else{
    modSummary<-rbind(modSummary,tmp)
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
  
  cat("Finished model:",i,"\n\n")
  
}

colnames(lassoCoeffs) <- varsToUse
rownames(lassoCoeffs) <- respVars
lassoCoeffs <- round(lassoCoeffs,3)
lassoCoeffs <- rbind(lassoCoeffs,Count=apply(lassoCoeffs,2,FUN = function(x) sum(x!=0)))

print(sort(apply(lassoCoeffs,2,FUN = function(x) sum(x!=0)), decreasing = TRUE))




write.csv(modSummary, "./OUT/modelsSummaries_v4z.csv")
write.csv(lassoCoeffs, "./OUT/modelSelectionCoeffsLasso_v4z.csv")

write.csv(vifRes, "./OUT/vifRes_v4z.csv")
write.csv(results, "./OUT/results_glm_v4z.csv", row.names = FALSE)


