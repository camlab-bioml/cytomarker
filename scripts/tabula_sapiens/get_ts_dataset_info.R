library(zellkonverter)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(plyr)
library(stringr)


ts_dir <- "/home/matt/cytomarker/tabula_sapiens/datasets/processed"
datasets <- list.files(ts_dir, ".rds", full.names = T)

dest_dir <- "/home/matt/cytomarker/tabula_sapiens/datasets/processed"

info_frame <- data.frame(tissue = character(), num_cells = numeric(),
                   num_genes = numeric(),
                   preview = character())
  
for (elem in datasets) {
  gc()
  tissue <- str_split_fixed(elem, "processed/|.rds", 3)[,2]
  tissue <- gsub(" ", "", tissue, fixed = TRUE)
  data <- readRDS(elem)
  print(tissue)

  compartments <- table(colData(data)$compartment) |> as.data.frame() |> mutate(pasted = paste(Var1, ": ",
                                                                                               Freq, " cells",
                                                                                               sep = ""))
  
  # len_possible_cats <- length(unique(colnames(colData(data))))
  # num_limit <- ifelse(len_possible_cats <= 3, len_possible_cats, 3)
  # cats_to_show <- unique(colnames(colData(data)))[1:num_limit]
  # others_addition <- ifelse(len_possible_cats <= 3, "", "and others")
  # var_group <- ifelse(len_possible_cats < 2, "metadata variable", "metadata variables")
  # info <- paste(len_possible_cats,
  #             var_group, "in selected dataset, including",
  #             paste(cats_to_show,
  #                   collapse = ", "),
  #             others_addition,
  #             sep = " ")
  
  info <- paste(length(unique(colData(data)$cell_ontology_class)), " subtypes in cell_ontology_class", 
                " <br> ", paste(compartments$pasted, collapse = " <br> "))
  
  
  to_add <- data.frame(tissue = tissue, num_cells = dim(data)[2],
                       num_genes = dim(data)[1],
                       preview = info)
  
  info_frame <- rbind(info_frame, to_add)
  
}

write.table(info_frame, file.path("inst", "ts_datasets.tsv"), 
                                row.names = F, sep = "\t", quote = F)

