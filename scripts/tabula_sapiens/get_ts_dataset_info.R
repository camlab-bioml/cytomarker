library(zellkonverter)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(plyr)
library(stringr)


ts_dir <- "/home/matt/cytosel/tabula_sapiens/datasets/processed"
datasets <- list.files(ts_dir, ".rds", full.names = T)

dest_dir <- "/home/matt/cytosel/tabula_sapiens/datasets/processed"

for (elem in datasets) {
  tissue <- str_split_fixed(elem, "TS_|.rds", 3)[,2]
  tissue <- gsub(" ", "", tissue, fixed = TRUE)
  data <- readRDS(elem)
  print(unique(data$cell_ontology_class))
}
