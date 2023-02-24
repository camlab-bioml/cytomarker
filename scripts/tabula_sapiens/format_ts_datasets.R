
library(zellkonverter)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(plyr)
library(stringr)

use_centroid <- FALSE

ts_dir <- "/home/matt/cytosel/tabula_sapiens/datasets"
datasets <- list.files(ts_dir, ".h5ad", full.names = T)
# datasets <- datasets[!grepl("Gland", datasets)]

dest_dir <- "/home/matt/cytosel/tabula_sapiens/datasets/processed"

for (data_point in datasets) {
  tissue <- str_split_fixed(data_point, "TS_|.h5ad", 3)[,2]
  tissue <- gsub(" ", "", tissue, fixed = TRUE)
  if (isFALSE(file.exists(file.path(dest_dir, paste(tissue, ".rds", sep = ""))))) {
  gc()
  print(data_point)
  data <- readH5AD(data_point)
  # data <- data[,data$n_counts_UMIs > 500 & data$n_genes > 500]
  
  umap <- reducedDim(data, "X_umap") |> as.data.frame() |> `colnames<-`(c("UMAP_1", "UMAP_2"))
  
  merged <- merge(colData(data), umap, by  = "row.names")
  
  merged <- merged |> as.data.frame() |> dplyr::group_by(cell_ontology_class) |> 
    dplyr::mutate(mean_UMAP_1 = median(UMAP_1),
                  sd_UMAP_1 = sd(UMAP_1),
                  mean_UMAP_2 = median(UMAP_2),
                  sd_UMAP_2 = sd(UMAP_2)) |> ungroup()
  
  
  types <- as.data.frame(table(data$cell_ontology_class)) |> filter(Freq > 20)
  
  types <- types |> mutate(subset_val = ifelse(Freq > 20 & Freq < 100, Freq, 
                               min(round_any(3000/length(types$Var1), 10), 
                                   min(types |> filter(Freq > 100) |> select(Freq))))) |>
    mutate(perc = 100*Freq/sum(Freq), frac_original = 100*(subset_val/Freq))
  
  # IMP: subset from all factors in the metadata column with the same number
  names <- c()
  for (elem in types$Var1) {
    new_data <- subset(merged, cell_ontology_class == elem)
    to_subset <- subset(types, Var1 == elem)$subset_val
    percent <- subset(types, Var1 == elem)$perc
    fraction <- subset(types, Var1 == elem)$frac_original
    
    if (fraction < 30 & isTRUE(use_centroid)) {
      
      deviation <- 1
      
      subset_umap <- new_data |> as.data.frame() |> filter((abs(UMAP_1) <= abs(mean_UMAP_1* + (deviation*sd_UMAP_1)) &
                                                              abs(UMAP_1) >= abs(mean_UMAP_1 + (1/deviation*sd_UMAP_1))) &
                                                             (abs(UMAP_2) <= abs(mean_UMAP_2 + (deviation*sd_UMAP_2)) &
                                                                abs(UMAP_2) >= abs(mean_UMAP_2 + (1/deviation*sd_UMAP_2))))
      
      # start in the middle of the cluster for each cell type and expand the search perimeter in
      # a breadth first sear manner until the proper number of cells is acquired
      while (to_subset > nrow(subset_umap)) {
        
        deviation <- 1.05*deviation
        subset_umap <- new_data |> as.data.frame() |> filter((abs(UMAP_1) <= abs(mean_UMAP_1* + (deviation*sd_UMAP_1)) &
                                                                abs(UMAP_1) >= abs(mean_UMAP_1 + (1/deviation*sd_UMAP_1))) &
                                                               (abs(UMAP_2) <= abs(mean_UMAP_2 + (deviation*sd_UMAP_2)) &
                                                                  abs(UMAP_2) >= abs(mean_UMAP_2 + (1/deviation*sd_UMAP_2))))
        
        print(tissue)
        print(elem)
        print(deviation)
        print(nrow(subset_umap))
        
      }
      
      new_data <- subset_umap
      
    }
    
    to_subset <- min(subset(types, Var1 == elem)$subset_val, nrow(new_data))
    new_data <- new_data |> sample_n(to_subset)
    names <- c(names, new_data$Row.names)
    
  }
  
  subset <- data[, colnames(data) %in% names]
  
  sce <- SingleCellExperiment(assays = list(logcounts = assay(subset)), 
                  reducedDims = SimpleList(UMAP = reducedDim(subset, "X_umap")), colData = colData(subset))
  sce$cell_ontology_class <- factor(sce$cell_ontology_class, levels = types$Var1)
  saveRDS(sce, file.path(dest_dir, paste(tissue, ".rds", sep = "")))
      }
  
}





