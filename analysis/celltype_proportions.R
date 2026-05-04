# ============================================================
# Calculate broad cell group proportions from a Seurat object
# ============================================================
# This script reads in a Seurat object, maps detailed cell-type
# annotations to broader biological groups, calculates the proportion
# of cells in each group, and saves the results as a CSV file.
#
# Usage:
# Rscript cell_group_proportions.R <OBJ> <dataset_name> <assay_name> <OUTPUT_DIR>
#
# Arguments:
# OBJ          Path to the input Seurat object (.rds)
# dataset_name Name of the dataset being analyzed
# assay_name   Assay type, e.g. "SC" or "SN"
# OUTPUT_DIR   Directory where the output CSV will be saved
# ============================================================


# -----------------------------
# command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
OBJ <- args[1]
dataset_name <- args[2]
assay_name <- args[3]
OUTPUT_DIR <- args[4]

# -----------------------------
# Load required packages
# -----------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(stringr)
library(purrr)

# -----------------------------
# Load Seurat object
# -----------------------------
obj <- readRDS(OBJ)

# ------------------------------------------------------------
# Define lookup table for grouping detailed cell-type labels
# into broader biological categories
# ------------------------------------------------------------
# The input Seurat metadata should contain a column called
# `combined_classification`, which stores the cell-type annotation
# for each cell/nucleus.
#
# Cell-type labels are cleaned using lowercase text and trimmed
# whitespace before matching to make the join more robust.
# ------------------------------------------------------------
celltype_lookup <- tribble(
  ~celltype, ~group,
  
  "Epithelial cells", "Epithelial",
  "Cholangiocytes", "Epithelial",
  "Neuroepithelial", "Epithelial",
  "Hepatocytes", "Epithelial",
  "Proximal tubule cells", "Epithelial",
  "Connecting tubule cells", "Epithelial",
  "Distal tubule cells", "Epithelial",
  "Pulmonary alveolar type I cells", "Epithelial",
  "Pulmonary alveolar type II cells", "Epithelial",
  "Clara cells", "Epithelial",
  #"Mesothelial cells", "Epithelial",
  "Ciliated cells", "Epithelial",
  "Secretory cells", "Epithelial",
  "Liver progenitor cell", "Epithelial",
  "Principal cells (Collecting duct system)", "Epithelial",
  "Ureteric Bud cells", "Epithelial",
  "α-intercalated cells (Collecting duct system)", "Epithelial",
  "β-intercalated cells (Collecting duct system)", "Epithelial",
  "aLOH", "Epithelial",
  "aLOH", "Epithelial",
  "PT", "Epithelial",
  "DCT", "Epithelial",
  "CD_IC_A", "Epithelial",
  "CD_PC", "Epithelial",
  "dLOH", "Epithelial",
  "CD_IC_B", "Epithelial",
  "CD_IC", "Epithelial",
  "Podo", "Epithelial",
  "Secretory cell", "Epithelial",
  "Loop of Henle cells", "Epithelial",
  "Inner Medullary cells (Collecting duct system)", "Epithelial",
  
  "Mesenchymal cells", "Mesenchymal",
  "Hepatic stellate cells", "Mesenchymal",
  "Mesangial cells", "Mesenchymal",
  "Stromal cells", "Mesenchymal",
  #"Adipocytes", "Mesenchymal",
  "Fibroblasts","Mesenchymal",
  
  "Mast cells", "Myeloid",
  "Macrophages", "Myeloid",
  "Macrophage", "Myeloid",
  "Microglia", "Myeloid",
  "Monocytes", "Myeloid",
  "Alveolar macrophages", "Myeloid",
  
  "T cells", "Lymphoid",
  "B cells", "Lymphoid",
  "NK cells", "Lymphoid",
  "NK", "Lymphoid",
  "Plasma cells", "Lymphoid",
  "Immune system cells", "Lymphoid",
    
  "Endothelial cells", "Endothelial",
  "Endothelial cell", "Endothelial",
  
  "Oligodendrocytes", "Glial",
  "Astrocytes", "Glial",
  
  "Neurons", "Neuronal",
  
  "Erythrocytes", "Erythroid"
) %>%
  mutate(celltype_clean = str_to_lower(str_trim(celltype)))

# ------------------------------------------------------------
# Calculate cell group proportions
# ------------------------------------------------------------
# 1. Add dataset and assay labels.
# 2. Clean the cell-type labels from the Seurat metadata.
# 3. Join detailed cell types to broader groups.
# 4. Assign unmatched cell types to "Other / ungrouped".
# 5. Count the number of cells/nuclei in each group.
# 6. Convert counts to proportions within each dataset-assay pair.
# ------------------------------------------------------------
prop_tbl <- obj@meta.data %>%
  mutate(
    Dataset = dataset_name,
    Assay = assay_name,
    celltype = combined_classification,
    celltype_clean = str_to_lower(str_trim(combined_classification))
  ) %>%
  left_join(celltype_lookup, by = "celltype_clean") %>%
  mutate(group = ifelse(is.na(group), "Other / ungrouped", group)) %>%
  dplyr::count(Dataset, Assay, group, name = "n") %>%
  group_by(Dataset, Assay) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ------------------------------------------------------------
# Identify cell types that did not match the lookup table
# ------------------------------------------------------------
# This is useful for checking whether additional cell-type labels
# should be added to the lookup table.
# ------------------------------------------------------------
meta <- obj@meta.data

meta$celltype <- meta$combined_classification
meta$celltype_clean <- stringr::str_to_lower(stringr::str_trim(meta$celltype))

check_unmatched <- meta %>%
  dplyr::left_join(celltype_lookup, by = "celltype_clean") %>%
  dplyr::filter(is.na(group))

unique_unmatched <- sort(unique(check_unmatched$celltype.x))

print(unique_unmatched)

# -----------------------------
# Save output table
# -----------------------------
# Output file contains:
# Dataset, Assay, group, n, and prop
# -----------------------------
write_csv(prop_tbl, paste0(OUTPUT_DIR,dataset_name,"_",assay_name,"_group_proportions.csv"))
