source("~/R_scripts/Wilcox_test.R")

args <- commandArgs(trailingOnly=TRUE) 
SC <- args[1] # 
SN <- args[2] #  
DATASET <- args[3] #Dataset name
SPECIES <- args[4] #mouse or human

print(paste0("Processing....", DATASET))

sc <- readRDS(SC)
sn <- readRDS(SN)

pattern <- ifelse(
  SPECIES == "human",
  "XIST",
  "xist"
)
sc_xist <- sc[["RNA"]]$data[pattern, ]
sn_xist <- sn[["RNA"]]$data[pattern, ]
Wilcox_and_AUC(sc_xist,sn_xist)
output_file <- paste0(DATASET,"_sc_xist.rds")
saveRDS(sc_xist, output_file)
output_file <- paste0(DATASET,"_sn_xist.rds")
saveRDS(sn_xist, output_file)


