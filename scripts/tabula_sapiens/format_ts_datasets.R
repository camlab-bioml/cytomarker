
library(zellkonverter)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(plyr)
library(stringr)


ts_dir <- "/home/matt/cytosel/tabula_sapiens/datasets"
datasets <- list.files(ts_dir, ".h5ad", full.names = T)

dest_dir <- "/home/matt/cytosel/tabula_sapiens/datasets/processed"

for (elem in datasets) {
  tissue <- str_split_fixed(elem, "TS_|.h5ad", 3)[,2]
  tissue <- gsub(" ", "", tissue, fixed = TRUE)
  if (isFALSE(file.exists(file.path(dest_dir, paste(tissue, ".rds", sep = ""))))) {
  gc()
  print(elem)
  data <- readH5AD(elem)
  # data <- data[,data$n_counts_UMIs > 500 & data$n_genes > 500]
  
  
  types <- as.data.frame(table(data$cell_ontology_class)) |> filter(Freq > 20)
  
  types <- types |> mutate(subset_val = ifelse(Freq > 20 & Freq < 100, Freq, 
                               min(round_any(3000/length(types$Var1), 10), 
                                   min(types |> filter(Freq > 100) |> select(Freq)))))
  
  print(types)
  
  # IMP: subset from all factors in the metadata column with the same number
  names <- c()
  for (elem in types$Var1) {
    to_subset <- subset(types, Var1 == elem)$subset_val
    sub_frame <- subset(as.data.frame(colData(data)), cell_ontology_class == elem) |> sample_n(to_subset)
    names <- c(names, rownames(sub_frame))
    
  }
  
  subset <- data[, colnames(data) %in% names]
  
  sce <- SingleCellExperiment(assays = list(counts = assay(subset, "raw_counts")),
                              colData = DataFrame(as.data.frame(colData(subset)) |> select(-free_annotation)))
  sce$cell_ontology_class <- factor(sce$cell_ontology_class, levels = types$Var1)
  
  sce <- logNormCounts(sce)
  sce <- SingleCellExperiment(assays = list(logcounts = assay(sce, "logcounts"),
                                            counts = assay(sce, "counts")),
                              colData = colData(sce))
  sce <- runUMAP(sce, exprs_values = "logcounts")
  saveRDS(sce, file.path(dest_dir, paste(tissue, ".rds", sep = "")))
      }
  
}





