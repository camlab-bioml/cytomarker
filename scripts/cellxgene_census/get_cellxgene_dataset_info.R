library(zellkonverter)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(plyr)
library(stringr)


cxg_dir <- "/home/matt/cytomarker/cellxgene_census"
datasets <- list.files(cxg_dir, ".rds", full.names = T)

info_frame <- data.frame(tissue = character(), num_cells = numeric(),
                         num_genes = numeric(),
                         preview = character(),
                         category_interest = character())

for (elem in c(datasets)) {
  gc()
  identifier <- str_split_fixed(elem, "cellxgene_census/|.rds", 3)[,2]
  data <- readRDS(elem)
  
  num_cell_cats <- length(unique(colData(data)$cell_type))
  
  if (num_cell_cats > 1) {
    info <- paste(num_cell_cats, " subtypes in cell_type", sep="")
    
    to_add <- data.frame(tissue = identifier, num_cells = dim(data)[2],
                         num_genes = dim(data)[1],
                         preview = info,
                         category_interest = c("cell_type"))
    
    info_frame <- rbind(info_frame, to_add)
    print(identifier)
  }
  
}

write.table(info_frame, file.path("inst", "other_datasets.tsv"), 
            row.names = F, sep = "\t", quote = F)

