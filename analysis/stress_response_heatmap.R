# =======================
# Loading input arguments 
# =======================
# FILEPATH: directory containing per dataset CSV results 
# OUTPUTDIR: directory to save heatmap output 
args <- commandArgs(trailingOnly = TRUE)
FILEPATH <- args[1]
OUTPUTDIR <- args[2]

# ======================
# Load and combined data
# ======================
# Read CSV files in per dataset, generated using the stress_response_calculation.R script
files <- list.files(path=FILEPATH, pattern = "*_stress_gene_level_results.csv", full.names = TRUE)
# Combine the files into one dataframe 
combined <- files %>% 
  lapply(read_csv) %>%
  bind_rows()

# =======================
# Load required libraries 
# =======================
library(dplyr)
library(tidyr)
library(pheatmap)
library(stringr)
library(tibble)
library(grid)

# ============================
# Standardize cell type labels 
# ============================
# Ensure consistent naming across datasets  
combined <- combined %>%
  mutate(
    celltype = recode(celltype,
                            "Endothelial cell" = "Endothelial cells",
                            "Macrophage" = "Macrophages",
                            "B cell" = "B cells",
                            "T cell" = "T cells",
                            "Hepatic stellate cell" = "Hepatic stellate cells",
                            "Cholangiocyte" = "Cholangiocytes",
                            "Fibroblast" = "Fibroblasts",
                            "Mesothelial cell" = "Mesothelial cells",
                            "Mesangial cell" = "Mesangial cells",
                            "Connecting tubule cell" = "Connecting tubule cells",
                            "Pulmonary alveolar type I cell" = "Pulmonary alveolar type I cells",
                            "Pulmonary alveolar type II cell" = "Pulmonary alveolar type II cells",
                            .default = celltype
    )
  )

# ==========================
# Prepare plotting dataframe
# ==========================
# create a unique column identifier combining dataset + cell type 
plot_df <- combined %>%
  mutate(col_id = paste(dataset, celltype, sep = " "))

# Keep only selected cell types for plotting 
plot_df <- plot_df %>%
  filter(celltype %in% c("T cells", "Macrophages", "Endothelial cells"))

# ===============================
# Define ordering of rows/columns 
# ===============================
celltype_order <- c("T cells", "Macrophages", "Endothelial cells")
dataset_order <- c("Andrews", "Wen", "Denisenko", "Koenitzer")

plot_df <- plot_df %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order),
    dataset = factor(dataset, levels = dataset_order)
  ) %>%
  arrange(celltype, dataset)

# ================================
# Define stress response gene list 
# ================================
# Human gene symbols 
stress_response_genes_human <- c(
  "FOSB","FOS","JUN","JUNB","JUND","ATF3","EGR1","HSPA1A","HSPA1B",
  "HSP90AB1","HSPA8","HSPB1","IER3","IER2","BTG1","BTG2","DUSP1",
  "MAFF","NR4A1","SOCS3","GADD45B","CYR61","CSRNP1","NFKBIZ","PPP1R15A",
  "ZFP36","GADD45G","NFIL3","CLDN4","KRT8","TNFRSF12A","FOSL2","PQLC1",
  "KLF6","ADAMTS1","ERF","PIM3","GM5244","LMNA","DDIT3","COQ10B",
  "KLF11","GM24276","AL480526","2900052L18RIK"
)
# Mouse gene symbols (converted from human format)
stress_response_genes_mouse <- str_to_title(str_to_lower(stress_response_genes_human))
# Combined list (used earlier for filtering if needed)
stress_response_genes_all <- c(stress_response_genes_human, stress_response_genes_mouse)

plot_df <- plot_df %>%
  mutate(gene = factor(gene, levels = stress_response_genes_all))

# =====================================
# Standardize gene names across species 
# =====================================
# Convert all gene names to upper case to merge human + mouse rows 
plot_df <- plot_df %>%
  mutate(gene_standard = str_to_upper(gene))

# Create mapping table for displaying both human and mouse gene names 
gene_labels <- data.frame(
  gene_standard = stress_response_genes_human,
  human = stress_response_genes_human,
  mouse = str_to_title(str_to_lower(stress_response_genes_human))
) %>%
  mutate(label = paste0(human, "/", mouse))

# Join labels onto main dataframe 
plot_df <- plot_df %>%
  left_join(gene_labels, by = "gene_standard")

# =====================
# Create heatmap matrix
# =====================
# Convert long data frome to matrix (genes x dataset/celltype combinations)
mat_auc <- plot_df %>%
  dplyr::select(gene_standard, col_id, AUC) %>%
  pivot_wider(names_from = col_id, values_from = AUC) %>%
  column_to_rownames("gene_standard") %>%
  as.matrix()

# Replace row names with "human/mouse" gene labels 
rownames(mat_auc) <- gene_labels$label[
  match(rownames(mat_auc), gene_labels$gene_standard)
]

# =================================
# Order columns and define metadata
# =================================
col_order_df <- plot_df %>%
  distinct(col_id, celltype, dataset) %>%
  arrange(celltype, dataset) %>%
  mutate(
    species = case_when(
      dataset %in% c("Andrews", "Wen") ~ "Human",
      dataset %in% c("Denisenko", "Koenitzer") ~ "Mouse"
    )
  )

# Reorder matrix columns 
mat_auc <- mat_auc[, col_order_df$col_id, drop = FALSE]

# Column labels (only dataset names)
labels_col <- as.character(col_order_df$dataset)

# ====================================
# Define gaps between cell type groups
# ====================================
celltype_sizes <- col_order_df %>%
  count(celltype) %>%
  pull(n)

gaps_col <- cumsum(celltype_sizes)
gaps_col <- gaps_col[-length(gaps_col)]

# ===================================
# Create column annotations (species)
# ===================================
annotation_col <- data.frame(
  Species = factor(col_order_df$species, levels = c("Human", "Mouse"))
)
rownames(annotation_col) <- col_order_df$col_id

annotation_colors <- list(
  Species = c(
    Human = "#7fc97f",
    Mouse = "#beaed4"
  )
)

# ===============================
# Define colour scale for heatmap
# ===============================
breaks <- seq(0.35, 0.95, length.out = 101)

# ============
# Save heatmap
# ============
png(filename = paste0(OUTPUTDIR, "stress_response_heatmap.png"),
    width = 8, height = 12, units = "in", res = 300)

# Draw heatmap 
plot <- pheatmap(
  mat_auc,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  color = colorRampPalette(c("#2b6cb0", "#f7f7f7", "#d53e4f"))(100),
  breaks = breaks,
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 10,
  labels_col = labels_col,
  angle_col = 45,
  cellwidth = 30,
  cellheight = 18,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_names_col = TRUE
)

# ================================================
# Add custom cell-type labels above heatmap blocks 
# ================================================
grid::grid.text(
  label = c("T cells", "Macrophages", "Endothelial cells"),
  x = unit(c(0.16, 0.32, 0.50), "npc"), # Manually tuned sitation 
  y = unit(0.96, "npc"),
  gp = gpar(fontsize = 12, fontface = "bold")
)

# Close file 
dev.off()


