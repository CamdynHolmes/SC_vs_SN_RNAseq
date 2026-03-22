source("~/R_scripts/Wilcox_test.R")

args <- commandArgs(trailingOnly=TRUE) 
SC <- args[1] # 
SN <- args[2] #  
DATASET <- args[3] #Dataset name
SPECIES <- args[4] #mouse or human

sc <- readRDS(SC)
sn <- readRDS(SN)

pattern <- ifelse(
  SPECIES == "human",
  "MALAT1",
  "malat1"
)
sc_malat1 <- sc[["RNA"]]$data[pattern, ]
sn_malat1 <- sn[["RNA"]]$data[pattern, ]
Wilcox_and_AUC(sc_malat1,sn_malat1)
output_file <- paste0(DATASET,"_sc_malat1.rds")
saveRDS(sc_malat1, output_file)
output_file <- paste0(DATASET,"_sn_malat1.rds")
saveRDS(sn_malat1, output_file)

pattern <- ifelse(
  SPECIES == "human",
  "NEAT1",
  "neat1"
)
sc_neat1 <- sc[["RNA"]]$data[pattern, ]
sn_neat1 <- sn[["RNA"]]$data[pattern, ]
Wilcox_and_AUC(sc_neat1,sn_neat1)
output_file <- paste0(DATASET,"_sc_neat1.rds")
saveRDS(sc_neat1, output_file)
output_file <- paste0(DATASET,"_sn_neat1.rds")
saveRDS(sn_neat1, output_file)