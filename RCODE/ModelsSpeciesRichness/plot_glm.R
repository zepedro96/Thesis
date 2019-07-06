library(dplyr)
library(ggplot2)

class_tp <- read.csv("./OUTtoShare/class_type_v1.csv")
glm_dt <- read.csv("./OUT/results_glm_v4z.csv")

glm_res <- left_join(class_tp, glm_dt, by = "respVar")
#glm_res <- select(glm_res, -3)

glm_res <- glm_res[-c(12:14),]

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
ggsave("./OUTtoShare/CoxSnell-v3.jpg", width = 10, height = 18, units = "cm" ) 

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
ggsave("./OUTtoShare/Nagelkerke-v3.jpg", width = 10, height = 18, units = "cm" ) 

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
ggsave("./OUTtoShare/Effron-v3.jpg", width = 10, height = 18, units = "cm" ) 

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
ggsave("./OUTtoShare/OdT.LRT-v3.jpg", width = 10, height = 18, units = "cm" ) 


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
ggsave("./OUTtoShare/OdT.DeanB-v3.jpg", width = 10, height = 18, units = "cm" ) 


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
ggsave("./OUTtoShare/OdT.DeanB2-v3.jpg", width = 10, height = 18, units = "cm" ) 
