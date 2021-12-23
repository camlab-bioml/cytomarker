#' Split columns into usable (good) and unusable (bad) categories
#' to be passed on to compute_fm
good_col <- function(sce, column) {
  
  good_or_bad <- list(good = c(), 
                      bad = list(colname = c(), n = c()))
  
  for(col in column) {
    n_unique_elements <- length(unique(colData(sce)[[col]]))
    
    if(n_unique_elements == 1 || n_unique_elements > 100) {
      good_or_bad$bad$colname <- c(good_or_bad$bad$colname, col)
      good_or_bad$bad$n <- c(good_or_bad$bad$n, n_unique_elements)
      
      good_or_bad$bad <- list(colname = good_or_bad$bad$colname,
                              n = good_or_bad$bad$n)
    } else {
      good_or_bad$good <- c(good_or_bad$good, col)
    }
  }
  
  good_or_bad <- list(good = good_or_bad$good,
                      bad = good_or_bad$bad)
  
  good_or_bad
}

#' Compute the findMarkers outputs and store
#' 
#' Note: this used to be part of get_markers but needs
#' split off after colouring
compute_fm <- function(sce, columns, pref_assay = "logcounts") {

  fms <- lapply(columns, function(col) {
    
    test_type <- ifelse(pref_assay == "counts", "binom", "t")
    fm <- findMarkers(sce, colData(sce)[[col]], test.type = test_type, assay.type = pref_assay)
  })
  
  names(fms) <- columns
  
  fms
}


#' Compute markers
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom scran findMarkers
get_markers <- function(fms, panel_size) {
  
  columns <- names(fms)

  marker <- list(recommended_markers = c(), scratch_markers = c(), top_markers = c())

  for(col in columns) {
    fm <- fms[[col]]
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
  }
  marker <- list(recommended_markers = recommended,
                   scratch_markers = scratch,
                   top_markers = top)
  marker
}
       
#' Given a list of top markers and the fms, get the associated
#' cell types
get_associated_cell_types <- function(markers, fms) {
  fm <- fms[[1]] # For now, we're only doing this for the first
  
  all_markers <- unique(unlist(markers[c('recommended_markers', 'scratch_markers', 'top_markers')]))
  
  associated_cell_types <- sapply(all_markers, function(tm) {
    pvals <- sapply(fm, function(f) {
      tmp <- f[tm,] # get summary statistics for this gene
      
      ## return p value if it's a marker, or 1 otherwise
      ifelse(tmp$summary.logFC < 0, 1, tmp$FDR)  
    })
    names(pvals)[which.min(pvals)]
  })
  associated_cell_types
}

#' Controlled way of setting current markers
#' Return them in a safe way!
set_current_markers_safely <- function(markers, fms) {
  markers$associated_cell_types <- get_associated_cell_types(markers, fms)
  
  # print(markers)
  markers
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
#' @importFrom tibble tibble
#' @importFrom naivebayes naive_bayes
#' @importFrom dplyr bind_rows
train_nb <- function(x,y, cell_types) {
  
  flds <- caret::createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)
  
  metrics <- lapply(flds, function(test_idx) {
  
    fit <- naive_bayes(x[-test_idx,], y[-test_idx])
    
    p <- predict(fit, newdata = x[test_idx,])
    
    overall <- yardstick::bal_accuracy_vec(y[test_idx], p)
    scores <- sapply(cell_types, function(ct) {
      yardstick::ppv_vec(factor(y[test_idx] == ct, levels=c("TRUE", "FALSE")), 
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
  expression_mat <- Heatmap(mat,
                            col = viridis(100),
                            name="Expression")
  
  if(display == "Marker-marker correlation") {
    return(cor_mat)
  } else {
    return(expression_mat)
  }
  
  
}

#' Suggest a set of redundant genes to remove
#' 
#' @param cmat An correlation matrix of *the full expression*  calculated as
#' expression <- as.matrix(assay(sce, pref_assay)[markers$top_markers,])
#' cmat <- cor(t(expression))
#' @param n_genes Number of genes to suggest to remove
#' @importFrom caret findCorrelation
suggest_genes_to_remove <- function(cmat, n_genes=10) {
  rg <- c()
  for(i in seq_len(n_genes)) {
    lgl <- !(rownames(cmat) %in% rg)
    fcs <- findCorrelation(cmat[lgl, lgl], cutoff = 0.01)
    gene_to_remove <- colnames(cmat[lgl, lgl])[ fcs[1] ]
    rg <- c(rg, gene_to_remove)
  }
  rg

}


round3 <- function(x) format(round(x, 1), nsmall = 3)


#' Read the input single-cell RNA-seq dataset
#' from compressed RDS file. Accepts either:
#' (1) SingleCellExperiment
#' (2) Seurat object (converted to SingleCellExperiment)
#' 
#' TODO: accept anndata
#' 
read_input_scrnaseq <- function(obj) {
  
  if(is(obj, 'SingleCellExperiment')) {
    return(obj)
  } else if(is(obj, 'Seurat')) {
    ## Convert to SingleCellExperiment
    sce <- Seurat::as.SingleCellExperiment(obj)
    return(sce)
  } 
} 

#' Return the legend of the cell type colours
#' @importFrom ggplot2 theme_bw geom_histogram scale_fill_manual guides
get_legend <- function(palette) {
  ## Draw the ggplot
  df <- tibble(r = runif(length(palette)),
               `Cell type` = names(palette))
  
  # TODO: Maybe "Cell type" isn't the best name for this
  plt <- ggplot(df, aes(x = r, fill = `Cell type`)) + 
    geom_histogram() +
    scale_fill_manual(values = palette) +
    guides(fill = guide_legend(ncol = 8)) 
  
  legend <- cowplot::get_legend(plt)
  legend
}

#' Map gene symbols to antibody icons depending on
#' whether the gene is found in the database
map_gene_name_to_antibody_icon <- function(gene_id, df_antibody) {
  antibody_info <- get_antibody_info(gene_id, df_antibody)
  # TODO: fix this
  # antibody_info <- list()
  # antibody_info$status <- sample(c("NO_GENE_FOUND", "NO_ANTIBODY_FOUND", "ANTIBODY_FOUND"), 1)
  
  if(antibody_info$status == "ANTIBODY_FOUND") {
    return(icon("check"))
  } else if(antibody_info$status == "NO_ANTIBODY_FOUND") {
    return(icon("times"))
  } else {
    return(icon("question"))
  }
  
}

#' Get antibody info/status for a particular gene
get_antibody_info <- function(gene_id, df) {

  
  df <- df[df$Symbol == gene_id,]
  
  if(nrow(df) == 0) {
    return(
      list(
        status = "NO_GENE_FOUND",
        antibodies = c()
      )
    )  
  } 
  
  ## This doesn't actually make sense with current version,
  ## will almost always return ANTIBODY_FOUND
  # status <- ifelse(is.na(df$ab_name), "NO_ANTIBODY_FOUND", "ANTIBODY_FOUND")
  status <- "ANTIBODY_FOUND"

  antibodies <- df$id

  
  list(
    status = status,
    antibodies = antibodies
  )
  
}

#' Compute alternative markers
compute_alternatives <- function(gene_to_replace, sce, pref_assay, n_correlations) {
  x <- as.matrix(assay(sce, pref_assay))[,sample(ncol(sce), min(5000, ncol(sce)))]
  
  y <- x[gene_to_replace, ]
  
  yo <- x[rownames(x) != gene_to_replace,]
  
  correlations <- cor(t(yo), y)
  
  alternatives <- data.frame(Gene = rownames(yo), Correlation = correlations[,1])
  alternatives <- alternatives[!is.na(alternatives$Correlation),]
  alternatives <- alternatives[order(-(alternatives$Correlation)),]
  alternatives <- alternatives[1:n_correlations,]
  rownames(alternatives) <- NULL
  
  alternatives
  
}

