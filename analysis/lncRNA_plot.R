############################################################
# Title: SC vs SN NEAT1/MALAT1 Expression and AUC Plot
# Author: Camdyn Holmes
# Description:
# This script reads RDS files containing gene expression values
# for single-cell (SC) and single-nucleus (SN) datasets, combines
# them into one dataframe, and creates a two-panel figure:
#   1) SC vs SN normalized expression boxplots
#   2) Dataset-specific AUC values
############################################################

#-----------------------------------------------------------
# Read command-line argument
#-----------------------------------------------------------

# DATA_DIR should be the folder containing the expression RDS files.
# Example file names are expected to look like:
#   Andrews_sc_malat1.rds
#   Andrews_sn_malat1.rds
#   Andrews_sc_neat1.rds
#   Andrews_sn_neat1.rds
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- args[1]

#-----------------------------------------------------------
# Load required libraries
#-----------------------------------------------------------
library(ggplot2)  # for plotting
library(cowplot)  # for combining plots
library(dplyr)    # for data wrangling
library(stringr)  # for parsing dataset/assay names from filenames
library(tibble)   # for creating AUC table

#-----------------------------------------------------------
# Read expression data from RDS files
#-----------------------------------------------------------
# Get full paths to all files in the input directory
data_files <- list.files(
  path = DATA_DIR, 
  full.names = TRUE
)

# Read each RDS file into a list
all_data <- lapply(data_files, readRDS)

# Name each list element using the file name without ".rds"
# These names are later used to extract dataset and assay information.
names(all_data) <- sub("\\.rds$", "", basename(data_files))

#-----------------------------------------------------------
# Define plotting order and colours
#-----------------------------------------------------------
# Fixed dataset order for consistent plotting across figures
ds_order <- c("Andrews","Bakken","Denisenko","Koenitzer",
              "Korrapati","Selewa","Valk","Wen","Wu")

# Colours for assay type
sc_col <- "#7fc97f"
sn_col <- "#beaed4"

# Dataset-specific colours for the AUC barplot
ds_cols <- c(
  Andrews   = "#8dd3c7",
  Bakken    = "darkblue",
  Denisenko = "#ffffb3",
  Koenitzer = "#bebada",
  Korrapati = "#fb8072",
  Selewa    = "#80b1d3",
  Valk      = "#fdb462",
  Wen       = "#b3de69",
  Wu        = "#fccde5"
)

##-----------------------------------------------------------
# Combine expression values into one dataframe
#-----------------------------------------------------------
df_box <- bind_rows(lapply(names(all_data), function(nm) {
  
  # Split file name into components.
  # Assumes file names follow the format:
  #   Dataset_Assay_Gene
  # Example:
  #   Andrews_sc_neat1
  parts <- str_split(nm, "_", simplify = TRUE)
  
  data.frame(
    Dataset = parts[1],
    Condition = toupper(parts[2]), # convert sc/sn to SC/SN
    Value = as.numeric(all_data[[nm]])
  )
}))

# Keep only expected datasets and assay types, then set factor order
df_box <- df_box %>%
  filter(Dataset %in% ds_order, Condition %in% c("SC", "SN")) %>%
  mutate(
    Dataset = factor(Dataset, levels = ds_order),
    Condition = factor(Condition, levels = c("SC", "SN"))
  )

#-----------------------------------------------------------
# Define dataset-specific AUC values
#-----------------------------------------------------------

# AUC values summarize how well expression separates SC and SN
# distributions for the gene being plotted.
#
# Interpretation:
#   AUC = 0.5 suggests little/no separation
#   AUC > 0.5 suggests stronger separation between SC and SN

# MALAT1 AUC values, if plotting MALAT1 instead of NEAT1:
# auc_tbl <- tibble(
#   Dataset = ds_order,
#   AUC = c(
#     0.6459618,0.55,0.9223013 , 0.9505932, 0.7450102, 0.7909317, 0.6045707, 0.6929307, 0.9048841
#   )
# )

# NEAT1 AUC values
auc_tbl <- tibble(
  Dataset = ds_order,
  AUC = c(
    0.7764101,0.53, 0.9248652, 0.6127468, 0.6017386, 0.505704, 0.5184642, 0.752845, 0.8347936
  )
)

# Match AUC table dataset order to plotting order
auc_tbl$Dataset <- factor(auc_tbl$Dataset, levels = ds_order)

#-----------------------------------------------------------
# Create AUC barplot
#-----------------------------------------------------------

p_auc <- ggplot(auc_tbl, aes(x = Dataset, y = AUC, fill = Dataset)) +
  coord_cartesian(ylim = c(0.5, 1.0)) + # Restrict visible y-axis range to focus on AUC values above 0.5
  geom_col(width = 0.7) + # Plot one bar per dataset
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "black", linewidth = 0.8) + # Reference line indicating moderate separation
  scale_fill_manual(values = ds_cols) +
  labs(x = "Dataset", y = "AUC") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 35),
    axis.text.x = element_text(size = 35, angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank()
  )

#-----------------------------------------------------------
# Create SC vs SN expression boxplot
#-----------------------------------------------------------
p_box <- ggplot(df_box, aes(x = Dataset, y = Value, fill = Condition)) +
  
  # Plot SC and SN expression distributions side-by-side for each dataset
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    outlier.shape = NA
  ) +
  
  scale_fill_manual(values = c(SC = sc_col, SN = sn_col)) +
  labs(x = NULL, y = "Normalized Expression", fill = NULL) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 35),
    axis.text.y  = element_text(size = 27),
    legend.text  = element_text(size = 45),
    legend.title = element_text(size = 40)
  )

#-----------------------------------------------------------
# Combine expression and AUC plots
#-----------------------------------------------------------
final_plot <- plot_grid(
  p_box,
  p_auc,
  ncol = 1
)

final_plot

#-----------------------------------------------------------
# Save final figure
#-----------------------------------------------------------
ggsave(
  final_plot,
  filename = "~/Documents/SC_vs_SN/poster_figures/neat1_plot.png",
  height = 13,
  width = 15,
  dpi = 300
)
