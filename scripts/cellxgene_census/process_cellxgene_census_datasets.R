library(cellxgene.census)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(data.table)
library(zellkonverter)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(plyr)
library(stringr)
library(Rtsne)

##### Get the list of all of the human datasets in the census ######

census <- open_soma()

# census_datasets <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())

# human_gene_metadata <- as.data.frame(census$get("census_data")$get("homo_sapiens")$obs$read()$concat())
# unique_human_lists <- c(unique(human_gene_metadata$dataset_id))
# fwrite(list(unique_human_lists), "cellxgene_census_human_datasets.csv")

# human_datasets <- read.csv("cellxgene_census_human_datasets.csv", header = F)
# 
# cleaned_datasets <- census_datasets |> filter(! collection_name %in% c("Tabula Sapiens") &
#                                                 dataset_id %in% human_datasets$V1 &
#                                                 dataset_total_cell_count < 100000 &
#                                                 dataset_total_cell_count >= 1000) |>
#   distinct(dataset_id, .keep_all = T)
# 
# write.table(cleaned_datasets, file.path("inst", "cellxgene_census_human_datasets.tsv"), quote = F, row.names = F,
#           sep = "\t")


dest_dir <- "/home/matt/cytomarker/cellxgene_census"

cleaned_datasets <- read.csv(file.path("inst", "cellxgene_census_human_datasets.tsv"), sep = "\t")

existing <- list.files(dest_dir, ".rds", full.names = F)

for (dataset in cleaned_datasets$dataset_id) {
  
  full_name <- paste(dataset, ".rds", sep="")
  if (!full_name %in% existing) {
    gc()
    print(dataset)
    query_ids <- c(dataset)
    sce <- get_single_cell_experiment(census, 
                                      organism = "Homo sapiens",
                                      obs_value_filter = "dataset_id %in% query_ids")
    
    sce_clean <- sce[rowSums(assay(sce)) > 10,]
    
    
    types <- as.data.frame(table(sce$cell_type)) |> filter(Freq > 20)
    
    types <- types |> mutate(subset_val = ifelse(Freq > 20 & Freq < 100, Freq, 
                                                 min(round_any(3000/length(types$Var1), 10), 
                                                     min(types |> filter(Freq > 100) |> select(Freq))))) |>
      mutate(perc = 100*Freq/sum(Freq), frac_original = 100*(subset_val/Freq))
    
    names <- c()
    for (elem in types$Var1) {
      new_data <- subset(as.data.frame(colData(sce)), cell_type == elem)
      to_subset <- subset(types, Var1 == elem)$subset_val
      percent <- subset(types, Var1 == elem)$perc
      fraction <- subset(types, Var1 == elem)$frac_original
      
      to_subset <- min(subset(types, Var1 == elem)$subset_val, nrow(new_data))
      new_data <- new_data |> sample_n(to_subset)
      names <- c(names, rownames(new_data))
      
    }
    
    subset <- sce_clean[, colnames(sce_clean) %in% names]
    subset$cell_type <- factor(subset$cell_type, levels = types$Var1)
    counts <- assay(subset)
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    logcounts(subset) <- log2(t(t(counts)/size.factors) + 1)
    subset <- scater::runUMAP(subset, exprs_values = "logcounts")
    sce <- SingleCellExperiment(assays = list(logcounts = assay(subset, "logcounts")), 
                                reducedDims = SimpleList(UMAP = reducedDim(subset)), 
                                colData = colData(subset))
    sce$cell_type <- factor(sce$cell_type, levels = types$Var1)
    if (dim(sce)[2] >= 1000 ) {
      saveRDS(sce, file.path(dest_dir, paste(dataset, ".rds", sep = "")))
    }
  }
}


