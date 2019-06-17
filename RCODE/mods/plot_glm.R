library(dplyr)
library(ggplot2)

class_tp <- read.csv("./OUT/class_type_v1.csv")
glm_dt <- read.csv("./OUT/results_glm_v2.csv")

glm_res <- left_join(class_tp, glm_dt, by = "respVar")
glm_res <- select(glm_res, -3)

tmp <- gsub("^Hab_","",glm_res$respVar)
tmp <- gsub("^Feed_","",tmp)
tmp <-  gsub("^Nest_","",tmp)
tmp <-  gsub("^Size_","",tmp)
tmp <-  gsub("_sum$","",tmp)

glm_res[,"respVar"] <- tmp

glm_res$respVar <- factor(glm_res$respVar, levels = glm_res$respVar[order(-glm_res$Effron)])

# 
# respVar <- select(glm_res, 1)
# CoxSnell <- select(glm_res, 3)
# Nagelkerke <- select(glm_res, 4)
# Effron <- select(glm_res, 5)
# varsgvargroups <- select(glm_res, 2)


p <- ggplot(glm_res, aes(x=respVar, y=Effron)) + 
  geom_bar(stat = "identity", fill="#7fbf7f", color="black") + 
  facet_wrap (~vargroups, scales = "free_x", ncol = 1) +
  xlab("Response Variable") +
  ylab("Effron R2") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

p <- plot(p)


