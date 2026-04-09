############################################################
# Script: calculate_mt_rb_metrics.R
#
# Description:
# This script compares single-cell (SC) and single-nucleus (SN)
# RNA-seq datasets by computing AUC values from Wilcoxon tests
# for:
#   1) Mitochondrial gene percentage (percent.mt)
#   2) Ribosomal gene percentage (percent.rb)
#
# Inputs:
#   1. SC (path to Seurat object for single-cell data)
#   2. SN (path to Seurat object for single-nucleus data)
#   3. species ("human" or "mouse")
#   4. DATASET_NAME (string used for output file naming)
#
# Output:
#   A CSV file containing AUC values for mitochondrial and
#   ribosomal gene percentages.
############################################################

############################
# 1. Parse command-line arguments
############################
args <- commandArgs(trailingOnly = TRUE)
SC <- args[1]
SN <- args[2]
species <- args[3] 
DATASET_NAME <- args[4]

############################
# 2. Load required functions and libraries
############################
# Custom function that performs Wilcoxon test and returns AUC
source("~/Wilcox_test.R.txt")

library(Seurat)

############################
# 3. Load Seurat objects
############################
sc <- readRDS(SC) # Single-cell dataset
sn <- readRDS(SN) # Single-nucleus dataset

############################
# 4. Extract mitochondrial percentages
############################
# Assumes percent.mt was already computed during preprocessing
sc_pct_mt <- sc$percent.mt 
sn_pct_mt <- sn$percent.mt

############################
# 5. Define ribosomal gene pattern based on species
############################
# Human genes: RPS / RPL
# Mouse genes: Rps / Rpl
if (species == "human") {
  pattern = "^RP[SL]"
} else if (species == "mouse") {
  pattern = "^Rp[sl]"
} else {
  print("Please select 'mouse' or 'human'")
}

############################
# 6. Calculate ribosomal gene percentages
############################
# Adds percent.rb metadata to Seurat objects
sc[["percent.rb"]] <- PercentageFeatureSet(sc, pattern = pattern)
sn[["percent.rb"]] <- PercentageFeatureSet(sn, pattern = pattern)

# Extract ribosomal percentages
sc_pct_rb <- sc$percent.rb
sn_pct_rb <- sn$percent.rb

############################
# 7. Compute AUC values using Wilcoxon test
############################
# Function returns multiple values; extract AUC (second element)
mt_auc <- Wilcox_and_AUC(sc_pct_mt,sn_pct_mt)[2]
rb_auc <- Wilcox_and_AUC(sc_pct_rb, sn_pct_rb)[2]

############################
# 8. Store results in a data frame
############################
results <- data.frame(
  dataset = DATASET_NAME,
  species = species,
  mt_auc = mt_auc,
  rb_auc = rb_auc
)

############################
# 9. Save results to CSV file
############################
write.csv(results, paste0(DATASET_NAME, "_auc_results.csv"), row.names = FALSE)




