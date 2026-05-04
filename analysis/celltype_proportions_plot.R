############################################################
# Title: Cell Type Group Proportion Plot (SC vs SN)
# Author: Camdyn Holmes
# Description:
# This script reads in cell type group proportion CSV files
# from multiple datasets, combines them into a single dataframe,
# and generates a stacked bar plot comparing proportions
# between single-cell (SC) and single-nucleus (SN) assays.
############################################################

# Load required libraries
library(readr) # for reading CSV files 
library(dplyr) # for data manipulation 
library(purrr) # for functinal programming (map)
library(ggplot2) # for plotting 
library(scales) # for formatting axes
library(stringr) # for extracting dataset/assay info from file names

#-----------------------------------------------------------
# Define input directory containing proportion CSV files
#-----------------------------------------------------------
# folder containing all your saved proportion CSVs
prop_dir <- "~/Documents/SC_vs_SN/Results/celltype_proportions"

files <- list.files(
  path = prop_dir,
  pattern = "_group_proportions.csv$",
  full.names = TRUE
)

#-----------------------------------------------------------
# Read and combine all files into a single dataframe
#-----------------------------------------------------------
plot_df <- map_dfr(files, function(f) {
  read_csv(f) %>%
    mutate(
      # store filename for reference 
      file = basename(f),
      
      # Extract dataset name (everything before first underscore)
      Dataset = str_extract(file, "^[^_]+"),
      
      # Extract assay type (SC or SN) from filename
      Assay = str_extract(file, "(?<=_)(SC|SN)(?=_)")
    )
})

#-----------------------------------------------------------
# Optional: set factor levels for consistent plotting order
#-----------------------------------------------------------

# Ensure datasets appear in a fixed order across plots
plot_df$Dataset <- factor(
  plot_df$Dataset,
  levels = c("Andrews", "Bakken","Denisenko", "Koenitzer",
             "Korrapati", "Selewa","Valk", "Wen", "Wu")
)

# Ensure SC is always plotted before SN
plot_df$Assay <- factor(plot_df$Assay, levels = c("SC", "SN"))

#-----------------------------------------------------------
# Create stacked bar plot
#-----------------------------------------------------------
p <- ggplot(plot_df, aes(x = Assay, y = prop, fill = group)) +
  
  # Stacked bars representing proportion of each cell type group
  geom_col(width = 0.75, colour = "black", linewidth = 0.15) +
  
  # Create one panel per dataset
  facet_wrap(~ Dataset, nrow = 1) +
  
  # Format y-axis as percentages
  scale_y_continuous(labels = percent_format()) +
  
  # Axis and legend labels
  labs(
    x = NULL,
    y = "Proportion of cells/nuclei",
    fill = "Cell type group"
  ) +
  theme_classic(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold", angle = 45),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )

p

#----------------------
# Save plot to file
#----------------------
ggsave(
  filename = "~/Documents/SC_vs_SN/poster_figures/celltype_group_stacked_barplot.png",
  plot = p,
  width = 18,
  height = 6,
  dpi = 600
)
