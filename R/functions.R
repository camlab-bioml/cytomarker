
#' Compute markers
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom scran findMarkers
get_markers <- function(sce, columns, panel_size, pref_assay = "logcounts") {
    
    marker <- list(recommended_markers = c(), scratch_markers = c(), top_markers = c())
  
    for(col in columns) {
      n_unique_elements <- length(unique(colData(sce)[[col]]))
      
      if(n_unique_elements == 1) { # need to fix this
        ## TODO: turn this into UI dialog
        stop(paste("Only one level in column", col))
      } else if(n_unique_elements > 100) {
        ## TODO: turn this into UI dialog
        stop(paste("Column", col, "has more than 100 unique elements"))        
      } else {
        test_type <- ifelse(pref_assay == "counts", "binom", "t")
        fm <- findMarkers(sce, colData(sce)[[col]], test.type = test_type, assay.type = pref_assay)
        
        n <- length(fm)
        
        top_select <- round(2 * panel_size / n)
        # all_select <- round(1000 / n)
        
        recommended <- marker$recommended_markers
        top <- marker$top_markers
        scratch <- marker$scratch_markers
        
        for(i in seq_len(n)) {
          f <- fm[[i]]
          recommended <- c(top, rownames(f)[seq_len(top_select)])
          top <- c(top, rownames(f)[seq_len(top_select)])
        }
        
        recommended <- unique(recommended)
        scratch <- unique(scratch)
        top <- unique(top)
        
        marker <- list(recommended_markers = recommended,
                       scratch_markers = scratch,
                       top_markers = top)
      }
    }
    
    
    marker

}

#' Compute the UMAP 
#' 
#' @importFrom scater runUMAP
#' @importFrom tibble tibble
#' @importFrom SingleCellExperiment reducedDim
get_umap <- function(sce, columns, pref_assay = "logcounts") {
  
  sce <- runUMAP(sce, exprs_values = pref_assay)
  
  df <- tibble(
    UMAP1 = reducedDim(sce, 'UMAP')[,1],
    UMAP2 = reducedDim(sce, 'UMAP')[,2]
  )
  
  for(column in columns) {
    df[[column]] <- as.data.frame(colData(sce))[[column]]
  }

  df
}

get_scores <- function(sce, columns, mrkrs, max_cells = 5000, pref_assay = "logcounts") {
  
  max_cells <- min(ncol(sce), max_cells)
  sce_tr <- sce[mrkrs, sample(ncol(sce), max_cells, replace=FALSE)]
  
  scores <- list()
  
  for(column in columns) {
    scores[[column]] <- get_scores_one_column(sce_tr, column, mrkrs, max_cells, pref_assay)
  }
  
  scores
}

get_scores_one_column <- function(sce_tr, column, mrkrs, max_cells = 5000, pref_assay = "logcounts") {
  
  x <- t(as.matrix(assay(sce_tr, pref_assay)))
  y <- factor(colData(sce_tr)[[column]])
  
  cell_types <- sort(unique(colData(sce_tr)[[column]]))
  
  train_nb(x,y, cell_types)
  
}

#' Train a Naive bayes classifier
#' 
#' @import caret
#' @importFrom tibble tibble
#' @import yardstick
#' @importFrom naivebayes naive_bayes
#' @importFrom dplyr bind_rows
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

#' Make the expression and correlation heatmaps
#' 
#' @importFrom ComplexHeatmap Heatmap 
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom stats cor
#' @importFrom viridis viridis
create_heatmap <- function(sce, markers, column, display, normalization, pref_assay = "logcounts") {

  normalization <- match.arg(normalization, c("Expression", "z-score"))
  
  if(display == "Marker-marker correlation") {
    mat <- summarizeAssayByGroup(
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
    
    # top_annot <- HeatmapAnnotation()
    cor_mat <- Heatmap(cor(t(mat)),
                       name="Correlation")
    
    return(cor_mat)
    
  } else if(display == column) {
    mat <- summarizeAssayByGroup(
      sce,
      id = colData(sce)[[display]],
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
                              name="Expression")
    
    return(expression_mat)
    
  }
  
}

round3 <- function(x) format(round(x, 1), nsmall = 3)


#' Read the input single-cell RNA-seq dataset
#' from compressed RDS file. Accepts either:
#' (1) SingleCellExperiment
#' (2) Seurat object (converted to SingleCellExperiment)
#' 
#' TODO: accept anndata
#' 
#' @importFrom Seurat as.SingleCellExperiment
read_input_scrnaseq <- function(input_path) {
  
  obj <- readRDS(input_path)
  
  if(is(obj, 'SingleCellExperiment')) {
    return(obj)
  } else if(is(obj, 'Seurat')) {
    ## Convert to SingleCellExperiment
    sce <- as.SingleCellExperiment(obj)
    return(sce)
  } else {
    ## TODO: Throw error: must be SingleCellExperiment or Seurat
    
  }
} 