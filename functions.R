
library(scran)
library(SingleCellExperiment)
library(tibble)
library(caret)
library(yardstick)
library(naivebayes)
library(ComplexHeatmap)
library(viridis)
library(dplyr)

# sce <- readRDS("data/sce.rds")


get_markers <- function(sce, column, panel_size, pref_assay = "logcounts") {
    test_type <- ifelse(pref_assay == "counts", "binom", "t")
    fm <- findMarkers(sce, colData(sce)[[column]], test.type = test_type, assay.type = pref_assay)

    n <- length(fm)
    
    top_select <- round(2 * panel_size / n)
    all_select <- round(1000 / n)
    
    recommended_markers <- c()
    top_markers <- c()
    all_markers <- c()
    
    for(i in seq_len(n)) {
      f <- fm[[i]]
      top_markers <- c(top_markers, rownames(f)[seq_len(top_select)])
      all_markers <- c(all_markers, rownames(f)[seq_len(all_select)])
      recommended_markers <- c(top_markers, rownames(f)[seq_len(top_select)])
    }
    
    top_markers <- unique(top_markers)
    all_markers <- unique(all_markers)
    recommended_markers <- unique(recommended_markers)
    
    all_markers <- setdiff(all_markers, top_markers)
    
    list(
      recommended_markers = recommended_markers,
      top_markers = top_markers,
      all_markers = all_markers
    )

}

get_umap <- function(sce, column, pref_assay = "logcounts") {
  sce <- runUMAP(sce, exprs_values = pref_assay)
  df <- tibble(
    UMAP1 = reducedDim(sce, 'UMAP')[,1],
    UMAP2 = reducedDim(sce, 'UMAP')[,2],
    x = colData(sce)[[column]]
  )
  names(df) <- c("UMAP1", "UMAP2", column)
  df
}

get_scores <- function(sce, column, mrkrs, max_cells = 5000, pref_assay = "logcounts") {
  max_cells <- min(ncol(sce), max_cells)
  sce_tr <- sce[mrkrs, sample(ncol(sce), max_cells, replace=FALSE)]
  
  x <- t(assay(sce_tr, pref_assay))
  x <- as.matrix(x)
  y <- factor(colData(sce_tr)[[ column ]])
  
  cell_types <- sort(unique(colData(sce)[[column]]))
  
  train_nb(x,y, cell_types)
  
}

train_nb <- function(x,y, cell_types) {
  flds <- createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)
  
  metrics <- lapply(flds, function(test_idx) {
  
    fit <- naive_bayes(x[-test_idx,], y[-test_idx])
    
    p <- predict(fit, newdata = x[test_idx,])
    
    overall <- bal_accuracy_vec(y[test_idx], p)
    scores <- sapply(cell_types, function(ct) {
      ppv_vec(factor(y[test_idx] == ct, levels=c("TRUE", "FALSE")), 
              factor(p == ct, levels = c("TRUE", "FALSE")))
    })
    tibble(
      what = c("Overall", cell_types),
      score=c(overall, scores)
    )
  }) %>% bind_rows()

  metrics
}

create_heatmap <- function(sce, markers, column, normalization, pref_assay = "logcounts") {
  normalization <- match.arg(normalization, c("Expression", "z-score"))
  
  mat <- scuttle::summarizeAssayByGroup(
    sce,
    id = colData(sce)[[column]],
    subset.row = markers$top_markers,
    statistics = 'mean',
    assay.type = pref_assay
  )
  mat <- (assay(mat, 'mean'))
  
  legend <- "Mean\nexpression"
  if(normalization == "z-score") {
    mat <- t(scale(t(mat)))
    legend <- "z-score\nexpression"
  }
  
  expression_mat <- Heatmap(mat,
                            col = viridis(100),
                            name="legend")
  cor_mat <- Heatmap(cor(t(mat)),
                     name="Correlation")
  expression_mat + cor_mat  
}

round3 <- function(x) format(round(x, 1), nsmall = 3)
