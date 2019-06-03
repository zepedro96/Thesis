library(readxl)
library(tidyr)
library(writexl)
library(dplyr)
#merge plots and data bird

codes <- read_xlsx("./DOCS/Birds_FieldData_Vez_ZP_v1.xlsx")

bird_trt <- read_xlsx("./DOCS/bd_traits_v2.xlsx")

vezDF <- merge(codes, bird_trt, by="Species_name")


#gather & sum

rsp <- gather(vezDF, 2:112, key = "code", value=presence, 
              convert = FALSE,na.rm = TRUE, factor_key = FALSE) %>% 
  mutate(ID_PSU = substr(.data$code,1,3)) %>%
  select(-code) %>% 
  filter(presence==1) %>% 
  group_by(ID_PSU, Species_name) %>% 
  summarize_all(.funs = first) %>%
  arrange(ID_PSU, Species_name) %>% 
  select(-presence, -Species_name) %>% 
  group_by(ID_PSU) %>% 
  summarize_all(.funs = list(sum=sum))


write_xlsx(rsp, "./DOCS/rsp.xlsx")


apply(rsp, 2, FUN = function(x) length(unique(x)) ) %>% sort
