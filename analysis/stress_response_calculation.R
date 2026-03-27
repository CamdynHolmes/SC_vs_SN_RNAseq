# ============================================
# Script: stress_response_calculation.R
# Author: Camdyn Holmes
#
# Description:
# This script compares stress-response gene expression between
# single-cell (SC) and single-nucleus (SN) Seurat objects.
#
# For each cell type shared between the SC and SN objects, it:
# 1. extracts expression values for a predefined list of
#    stress-response genes,
# 2. compares SC vs SN expression for each gene using a
#    Wilcoxon test + AUC calculation,
# 3. calculates mean expression and log2 fold change,
# 4. writes the results to a CSV file.
#
# Inputs:
#   args[1] = path to SC Seurat object (.rds)
#   args[2] = path to SN Seurat object (.rds)
#   args[3] = species ("human" or "mouse")
#   args[4] = dataset name
#   args[5] = output directory
#
# Output:
#   A CSV file with one row per gene per cell type containing:
#   dataset, celltype, gene, mean_sc, mean_sn, log2FC, AUC, p_value
# ============================================

# =======================
# Loading input arguments 
# =======================
args <- commandArgs(trailingOnly = TRUE)
SC <- args[1]
SN <- args[2]
species <- args[3]
dataset_name <- args[4]
output_dir <- args[5]

# ==============
# Load libraries 
# ==============
library(Seurat)
library(stringr)

# ============================
# Load custom functions
# ============================
# Wilcox_and_AUC(): performs Wilcoxon rank-sum test and computes AUC
# Make sure this path is correct on your system OR replace with a relative path
source("~/Downloads/Wilcox_test.R")

# ===================
# Read Seurat Objects
# ===================
sc <- readRDS(SC)
sn <- readRDS(SN)

# ============================
# Define stress response genes
# ============================
# Start with human-format gene list 
# If the species is mouse, convert to mouse-style gene names 
if (species == "human"){
  
  stress_response_genes <- c(
    "FOSB","FOS","JUN","JUNB","JUND","ATF3","EGR1","HSPA1A","HSPA1B",
    "HSP90AB1","HSPA8","HSPB1","IER3","IER2","BTG1","BTG2","DUSP1",
    "MAFF","NR4A1","SOCS3","GADD45B","CYR61","CSRNP1","NFKBIZ","PPP1R15A",
    "ZFP36","GADD45G","NFIL3","CLDN4","KRT8","TNFRSF12A","FOSL2","PQLC1",
    "KLF6","ADAMTS1","ERF","PIM3","GM5244","LMNA","DDIT3","COQ10B",
    "KLF11","GM24276","AL480526","2900052L18RIK"
  )
  
} else if (species == "mouse") {
  stress_response_genes <- c(
    "FOSB","FOS","JUN","JUNB","JUND","ATF3","EGR1","HSPA1A","HSPA1B",
    "HSP90AB1","HSPA8","HSPB1","IER3","IER2","BTG1","BTG2","DUSP1",
    "MAFF","NR4A1","SOCS3","GADD45B","CYR61","CSRNP1","NFKBIZ","PPP1R15A",
    "ZFP36","GADD45G","NFIL3","CLDN4","KRT8","TNFRSF12A","FOSL2","PQLC1",
    "KLF6","ADAMTS1","ERF","PIM3","GM5244","LMNA","DDIT3","COQ10B",
    "KLF11","GM24276","AL480526","2900052L18RIK"
  )
  #Convert to mouse format
  stress_response_genes <- str_to_title(str_to_lower(stress_response_genes))
} else{
  message(paste(species, "not included"))
}

# =========================================
# Keep genes only present in both SC and SN 
# =========================================
# First, keep stress response genes present in the SC object 
# Then further restrict to genes also present in the SN object 
common_stress_genes <- intersect(stress_response_genes, rownames(LayerData(sc, assay = "RNA", layer = "data")))
common_stress_genes <- intersect(common_stress_genes, rownames(LayerData(sn,assay = "RNA", layer = "data")))

# report stress-response genes that are missing from one or both objects 
missing_genes <- setdiff(stress_response_genes, common_stress_genes)
cat("Missing stress genes:\n")
print(missing_genes)
cat("Number missing:", length(missing_genes), "\n")

# ==================================
# Extract stress expression matrices 
# ==================================
# These matrices only contain genes of interest, 
# but still include all cells in the object
sc_stress_genes <- LayerData(sc, assay = "RNA",layer = "data")[common_stress_genes, , drop = FALSE]
sn_stress_genes <- LayerData(sn,assay = "RNA",layer = "data")[common_stress_genes, , drop = FALSE]

# ======================
# Find shared cell types 
# ======================
# Only compare cell types present in both SC and SN objects 
common_cts   <- intersect(sc$combined_classification, sn$combined_classification)

# Initialize list to store results for each cell types 
results_list <- list()

# Minimum number of cells required in each group
# This prevents the loop from failing when a cell type has low number of cells 
min_cells <- 3 

# ===========================
# Loop over shared cell types 
# ===========================
for (celltype in common_cts) {
  cat("Processing:", celltype, "\n")
  
  # subset cells by cell type
  sc_cells <- subset(sc, subset = combined_classification == celltype)
  sn_cells <- subset(sn, subset = combined_classification == celltype)
  
  cat("SC cells:", ncol(sc_cells), "\n")
  cat("SN cells:", ncol(sn_cells), "\n")
  
  # Skip this cell type if either SC or SN has too few cells 
  if (ncol(sc_cells) < min_cells || ncol(sn_cells) < min_cells) {
    cat("Skipping", celltype, "- too few cells in one group\n")
    next
  }
  
  # Extract normalized expression matrices for this cell type 
  sc_mat <- LayerData(sc_cells, layer = "data")
  sn_mat <- LayerData(sn_cells, layer = "data")
  
  # Store per-gene results for this cell type 
  gene_results <- list()
  
  # ===============================
  # Loop over stress-response genes 
  # ===============================
  for (gene in common_stress_genes) {
    
    # Skip if the gene is missing from either matrix 
    if (!(gene %in% rownames(sc_mat)) || !(gene %in% rownames(sn_mat))) {
      next
    }
    
    # Extract per-cell expression values for this gene   
    sc_vals <- as.numeric(sc_mat[gene, ])
    sn_vals <- as.numeric(sn_mat[gene, ])
    
    # Run Wilcoxon test and calculate AUC 
    wilcox_result <- Wilcox_and_AUC(sc_vals, sn_vals)
    
    # Store summary statistics for this gene in this cell type 
    gene_df <- data.frame(
      dataset = dataset_name,
      celltype = celltype,
      gene = gene,
      mean_sc = mean(sc_vals),
      mean_sn = mean(sn_vals),
      log2FC = log2((mean(sc_vals) + 1e-6) / (mean(sn_vals) + 1e-6)),
      AUC = wilcox_result[2],
      p_value = wilcox_result[4]
    )
    
    gene_results[[gene]] <- gene_df
  }
  
  # Combined all genes for this cell type into one dataframe 
  results_list[[celltype]] <- do.call(rbind, gene_results)
}

# =============================
# Combine all cell type results 
# =============================
results_df <- do.call(rbind, results_list)

# Print first few rows for quick checking 
print(head(results_df))

# ====================
# Write results to CSV 
# ====================
write.csv(
  results_df,
  paste0(output_dir, "/", dataset_name, "_stress_gene_level_results.csv"),
  row.names = FALSE
)







