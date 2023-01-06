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

info_frame <- data.frame(tissue = character(), num_cells = numeric(),
                   num_genes = numeric(),
                   preview = character())
  
for (elem in datasets) {
  gc()
  tissue <- str_split_fixed(elem, "processed/|.rds", 3)[,2]
  tissue <- gsub(" ", "", tissue, fixed = TRUE)
  data <- readRDS(elem)
  print(tissue)
  len_possible_cats <- length(unique(data[["cell_ontology_class"]]))
  num_limit <- ifelse(len_possible_cats <= 3, len_possible_cats, 3)
  cats_to_show <- unique(data[["cell_ontology_class"]])[1:num_limit]
  others_addition <- ifelse(len_possible_cats <= 3, "", "and others")
  var_group <- ifelse(len_possible_cats < 2, "grouping", "groupings")
  info <- paste(len_possible_cats,
              var_group, "in selected category, including",
              paste(cats_to_show,
                    collapse = ", "),
              others_addition,
              sep = " ")
  
  to_add <- data.frame(tissue = tissue, num_cells = dim(data)[2],
                       num_genes = dim(data)[1],
                       preview = info)
  
  info_frame <- rbind(info_frame, to_add)
  
}

write.csv(info_frame, "ts_datasets.csv", row.names = F, quote = F)

