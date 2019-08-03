

library(dplyr)
library(ggplot2)

# Changed the groups of variables because of modifications in names
class_tp <- read.csv("./OUTtoShare/class_type_v2.csv")

# Load results
glm_dt <- read.csv("./OUT/results_glm_elasticnet_v5.csv")

# Join tables
glm_res <- left_join(class_tp, glm_dt, by = "respVar")
#glm_res <- select(glm_res, -3)

# Remove variables related to feeding/omnivorous group
# Kept the one with best results
glm_res <- glm_res[-c(13:14),]

# Remove name prefixes
tmp <- gsub("^Hab_","",glm_res$respVar)
tmp <- gsub("^Feed_","",tmp)
tmp <-  gsub("^Nest_","",tmp)
tmp <-  gsub("^Size_","",tmp)
tmp <-  gsub("_sum$","",tmp)
glm_res[,"respVar"] <- tmp

# 
# respVar <- select(glm_res, 1)
# CoxSnell <- select(glm_res, 3)
# Nagelkerke <- select(glm_res, 4)
# Effron <- select(glm_res, 5)
# varsgvargroups <- select(glm_res, 2)

#-------------------------------------------------
#CoxSnell

glm_res$respVar <- factor(glm_res$respVar, levels = glm_res$respVar[order(-glm_res$CoxSnell)])


p <- ggplot(glm_res, aes(x=respVar, y=CoxSnell)) + 
  geom_bar(stat = "identity", fill="#7fbf7f", color="black") + 
  facet_wrap (~vargroups, scales = "free_x", ncol = 1) +
  xlab("Response Variable") +
  ylab("CoxSnell R2") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
  
p <- plot(p)
ggsave("./OUTtoShare/CoxSnell-v4.jpg", width = 10, height = 18, units = "cm" ) 

#-------------------------------------------------
#Nagelkerke

glm_res$respVar <- factor(glm_res$respVar, levels = glm_res$respVar[order(-glm_res$Nagelkerke)])


p <- ggplot(glm_res, aes(x=respVar, y=Nagelkerke)) + 
  geom_bar(stat = "identity", fill="#7fbf7f", color="black") + 
  facet_wrap (~vargroups, scales = "free_x", ncol = 1) +
  xlab("Response Variable") +
  ylab("Nagelkerke R2") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

p <- plot(p)
ggsave("./OUTtoShare/Nagelkerke-v4.jpg", width = 10, height = 18, units = "cm" ) 

#-------------------------------------------------
#Effron

glm_res$respVar <- factor(glm_res$respVar, levels = glm_res$respVar[order(-glm_res$Effron)])


p <- ggplot(glm_res, aes(x=respVar, y=Effron)) + 
  geom_bar(stat = "identity", fill="#7fbf7f", color="black") + 
  facet_wrap (~vargroups, scales = "free_x", ncol = 1) +
  xlab("Response Variable") +
  ylab("Effron R2") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

p <- plot(p)
ggsave("./OUTtoShare/Effron-v4.jpg", width = 10, height = 18, units = "cm" ) 

#-------------------------------------------------
#OdT.LRT

glm_res$respVar <- factor(glm_res$respVar, levels = glm_res$respVar[order(-glm_res$OdT.LRT)])


p <- ggplot(glm_res, aes(x=respVar, y=OdT.LRT)) + 
  geom_bar(stat = "identity", fill="#7fbf7f", color="black") + 
  facet_wrap (~vargroups, scales = "free_x", ncol = 1) +
  xlab("Response Variable") +
  ylab("OdT.LRT") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

p <- plot(p)
ggsave("./OUTtoShare/OdT.LRT-v4.jpg", width = 10, height = 18, units = "cm" ) 


#-------------------------------------------------
#OdT.DeanB

glm_res$respVar <- factor(glm_res$respVar, levels = glm_res$respVar[order(-glm_res$OdT.DeanB)])


p <- ggplot(glm_res, aes(x=respVar, y=OdT.DeanB)) + 
  geom_bar(stat = "identity", fill="#7fbf7f", color="black") + 
  facet_wrap (~vargroups, scales = "free_x", ncol = 1) +
  xlab("Response Variable") +
  ylab("OdT.DeanB") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

p <- plot(p)
ggsave("./OUTtoShare/OdT.DeanB-v4.jpg", width = 10, height = 18, units = "cm" ) 


#-------------------------------------------------
#OdT.DeanB2

glm_res$respVar <- factor(glm_res$respVar, levels = glm_res$respVar[order(-glm_res$OdT.DeanB2)])


p <- ggplot(glm_res, aes(x=respVar, y=OdT.DeanB2)) + 
  geom_bar(stat = "identity", fill="#7fbf7f", color="black") + 
  facet_wrap (~vargroups, scales = "free_x", ncol = 1) +
  xlab("Response Variable") +
  ylab("OdT.DeanB") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

p <- plot(p)
ggsave("./OUTtoShare/OdT.DeanB2-v4.jpg", width = 10, height = 18, units = "cm" )


#-------------------------------------------------
# Variable selection frequency
modSelFreq <- read.csv("./OUT/varSelectionFrequency_elasticnet_v5.csv")

modSelFreq$varName <- factor(modSelFreq$varName, 
                             levels = modSelFreq$varName[order(modSelFreq$freq, decreasing = TRUE)])

p <- ggplot(modSelFreq, aes(x=varName, y=relfreq)) + 
  geom_bar(stat = "identity", fill="#7fbf7f", color="black") + 
  xlab("Predictive Variable") +
  ylab("Relative frequency of selection\nin elasticnet models (%)") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

plot(p)
ggsave("./OUTtoShare/varSelectionFreq-ElasticNet-v4.jpg", plot = p, width = 14, height = 10, units = "cm" )


#-------------------------------------------------
# Average coeffs
lassoCoeffs <- read.csv("./OUT/ModelSelection_CoeffsElasticnet_v5.csv")

medCoeff <- data.frame(varName = colnames(lassoCoeffs)[-1], 
                       med=apply(lassoCoeffs[,-1], 2, median))


medCoeff$varName <- factor(medCoeff$varName, 
                             levels = medCoeff$varName[order(medCoeff$med, decreasing = TRUE)])

write.csv(medCoeff, "./OUT/MedianElasticnetCoeffs_v5.csv", row.names = FALSE)


p <- ggplot(medCoeff, aes(x=varName, y=med)) + 
  geom_bar(stat = "identity", fill="#7fbf7f", color="black") + 
  xlab("Predictive Variable") +
  ylab("Median coefficient value (elasticnet)") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

plot(p)
ggsave("./OUTtoShare/MedianCoeffValue-ElasticNet-v4.jpg", plot = p, width = 14, height = 10, units = "cm" )




