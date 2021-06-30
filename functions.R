
library(scran)
library(SingleCellExperiment)
library(tibble)

sce <- readRDS("data/sce.rds")


get_markers <- function(sce, column, panel_size) {

    
    fm <- findMarkers(sce, colData(sce)[[column]])

    
    n <- length(fm)
    
    top_select <- round(panel_size / n)
    all_select <- round(100 / n)
    
    top_markers <- c()
    all_markers <- c()
    
    for(i in seq_len(n)) {
      print(i)
      f <- fm[[i]]
      top_markers <- c(top_markers, rownames(f)[seq_len(top_select)])
      all_markers <- c(all_markers, rownames(f)[seq_len(all_select)])
    }
    
    top_markers <- unique(top_markers)
    all_markers <- unique(all_markers)
    
    all_markers <- setdiff(all_markers, top_markers)
    
    list(
      top_markers = top_markers,
      all_markers = all_markers
      
    )

}

get_umap <- function(sce, column) {
  sce <- runUMAP(sce)
  df <- tibble(
    UMAP1 = reducedDim(sce, 'UMAP')[,1],
    UMAP2 = reducedDim(sce, 'UMAP')[,2],
    x = colData(sce)[[column]]
  )
  names(df) <- c("UMAP1", "UMAP2", column)
  df
}




