
args <- commandArgs(trailingOnly = TRUE)
SC <- args[1]
SN <- args[2]
species <- args[3]
dataset_name <- args[4]
output_dir <- args[5]

library(Seurat)
library(stringr)

sc <- readRDS(SC)
sn <- readRDS(SN)

#Function to calculate the avg exp per gene per cell type 
avg_by_group <- function(seurat_obj,
                         group_col = "combined_classification", #the meta data column where your cell type annotation are 
                         assay = "RNA",
                         slot = "data") {
  # Set identity to grouping variable
  Idents(seurat_obj) <- seurat_obj[[group_col, drop = TRUE]]
  # Get average expression
  avg_mat <- AverageExpression(
    seurat_obj,
    assay = assay,
    slot = slot,
    return.seurat = FALSE
  )
  return(avg_mat[[assay]])
}

#compute per cell type averages 
sc_avg <- avg_by_group(sc)
sn_avg <- avg_by_group(sn)

#Keep only shared genes + shared cell types 
common_genes <- intersect(rownames(sc_avg), rownames(sn_avg))
common_cts   <- intersect(colnames(sc_avg), colnames(sn_avg))


sc_avg <- sc_avg[common_genes, common_cts, drop = FALSE]
sn_avg <- sn_avg[common_genes, common_cts, drop = FALSE]


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
  message(paste(species, "does not exist"))
}

# Keep only stress genes present in BOTH SC and SN
common_stress_genes <- intersect(stress_response_genes, rownames(sc_avg))
common_stress_genes <- intersect(common_stress_genes, rownames(sn_avg))

# print missing genes
missing_genes <- setdiff(stress_response_genes, common_stress_genes)
cat("Missing stress genes:\n")
print(missing_genes)
cat("Number missing:", length(missing_genes), "\n")

sc_stress_genes <- sc_avg[common_stress_genes, common_cts, drop = FALSE]
sn_stress_genes <- sn_avg[common_stress_genes, common_cts, drop = FALSE]

# Run Wilcoxon + AUC for each shared cell type
results_list <- vector("list", length(common_cts))
names(results_list) <- common_cts

for (celltype in common_cts) {
  celltype_avg_sc <- as.numeric(sc_stress_genes[, celltype])
  celltype_avg_sn <- as.numeric(sn_stress_genes[, celltype])
  
  res <- Wilcox_and_AUC(celltype_avg_sc, celltype_avg_sn)
  
  results_list[[celltype]] <- data.frame(
    dataset = dataset_name,
    celltype = celltype,
    n_genes = length(common_stress_genes),
    p_value = res[4],
    auc = res[2],
    stringsAsFactors = FALSE
  )
}

results_df <- do.call(rbind, results_list)

print(results_df)

write.csv(results_df, paste0(output_dir, dataset_name, "_stress_auc_results.csv"), row.names = FALSE)






