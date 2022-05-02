#' Split columns into usable (good) and unusable (bad) categories
#' to be passed on to compute_fm
#' 
#' @param sce A SingleCellExperiment object
#' @param columns A vector storing columns
good_col <- function(sce, columns) {
  
  good_or_bad <- list(good = c(), 
                      bad = list(colname = c(), n = c()))
  
  for(col in columns) {
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
#' 
#' @param sce A SingleCellExperiment object
#' @param columns A vector storing columns
#' @param pref_assay Assay loaded, default is 'logcounts'
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom scran findMarkers
compute_fm <- function(sce, columns, pref_assay, allowed_genes) {


  fms <- lapply(columns, function(col) {
    
    test_type <- ifelse(pref_assay == "counts", "binom", "t")
    fm <- findMarkers(sce, colData(sce)[[col]], test.type = test_type, assay.type = pref_assay)
    
    for(n in names(fm)) {
      fm[[n]] <- fm[[n]][rownames(fm[[n]]) %in% allowed_genes,]
    }
    fm
  })
  
  names(fms) <- columns
  
  fms
}


#' Compute markers
#' 
#' @param fms Stored findMarkers outputs
#' @param panel_size Targeted number of markers in selected panel
#' @param marker_strategy Strategy for picking markers, either heuristic 
#' cell type based or geneBasis
#' @param sce SingleCellExperiment object
#' @import geneBasisR
get_markers <- function(fms, panel_size, marker_strategy, sce) {
  
  columns <- names(fms)

  marker <- list(recommended_markers = c(), scratch_markers = c(), top_markers = c())
  
  if(marker_strategy == "geneBasis") {
    sce2 <- retain_informative_genes(sce[allowed_genes,])
    genes <- gene_search(sce2, n_genes=panel_size)
    marker <- list(
      recommended_markers = genes$gene,
      scratch_markers = c(),
      top_markers = genes$gene
    )
  } else {
    for(col in columns) {
      fm <- fms[[col]]
      n <- length(fm)
        
      top_select <- round(panel_size / n)
      
      recommended <- marker$recommended_markers
      top <- marker$top_markers
      scratch <- marker$scratch_markers
        
      for(i in seq_len(n)) {
        f <- fm[[i]]
        f <- f[!(rownames(f) %in% recommended),]
        
        ## Only keep markers that are over-expressed
        f <- f[f$summary.logFC > 0,]
        
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
  }

  marker
}
       
#' Given a list of top markers and the fms, get the associated
#' cell types
#' 
#' @param markers A list storing three lists of markers: 
#' recommended_markers, scratch_markers, and top_markers
#' @param fms Stored findMarkers outputs
get_associated_cell_types <- function(markers, fms) {
  fm <- fms[[1]] # For now, we're only doing this for the first
  
  all_markers <- unique(unlist(markers[c('recommended_markers', 'scratch_markers', 'top_markers')]))
  
  associated_cell_types <- sapply(all_markers, function(tm) {
    pvals <- sapply(fm, function(f) {
      tmp <- f[tm,] # get summary statistics for this gene
      
      ## return p value if it's a marker, or 1 otherwise
      # ifelse(tmp$summary.logFC < 0, 1, tmp$FDR)
      ## New: return summary logFC
      -tmp$summary.logFC
    })
    names(pvals)[which.min(pvals)]
  })
  associated_cell_types
}

#' Controlled way of setting current markers
#' Return them in a safe way!
#' 
#' @param markers A list storing three lists of markers: 
#' recommended_markers, scratch_markers, and top_markers
#' @param fms Stored findMarkers outputs
set_current_markers_safely <- function(markers, fms) {
  markers$associated_cell_types <- get_associated_cell_types(markers, fms)

  markers
}


#' Compute the UMAP 
#'
#' @param sce A SingleCellExperiment object
#' @param columns A vector storing columns
#' @param pref_assay Assay loaded, default is 'logcounts'
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

#' Compute the metric scores
#' 
#' @param sce A SingleCellExperiment object
#' @param columns A vector storing columns
#' @param mrkrs A list storing three lists of markers: 
#' recommended_markers, scratch_markers, and top_markers
#' @param max_cells Maximum number of cells
#' @param pref_assay Assay loaded, default is 'logcounts'
get_scores <- function(sce, columns, mrkrs, max_cells = 5000, pref_assay = "logcounts") {
  
  max_cells <- min(ncol(sce), max_cells)
  sce_tr <- sce[mrkrs, sample(ncol(sce), max_cells, replace=FALSE)]
  
  scores <- list()
  
  for(column in columns) {
    scores[[column]] <- get_scores_one_column(sce_tr, column, mrkrs, max_cells, pref_assay)
  }
  
  scores
}

#' Compute the metric scores for a specific column
#' 
#' @param sce_tr A SingleCellExperiment object
#' @param column A vector storing a column
#' @param mrkrs A list storing three lists of markers: 
#' recommended_markers, scratch_markers, and top_markers
#' @param max_cells Maximum number of cells
#' @param pref_assay Assay loaded, default is 'logcounts'
get_scores_one_column <- function(sce_tr, column, mrkrs, max_cells = 5000, pref_assay = "logcounts") {
  
  x <- t(as.matrix(assay(sce_tr, pref_assay)))
  y <- factor(colData(sce_tr)[[column]])
  
  cell_types <- sort(unique(colData(sce_tr)[[column]]))
  
  train_nb(x,y, cell_types)
}

#' Train a Naive bayes classifier
#' 
#' @param x A matrix calculated as `x <- t(as.matrix(assay(sce_tr)[[column]]))`
#' @param y A factor calculated as `y <- factor(colData(sce_tr)[[column]])`
#' @param cell_types A sorted vector with unique elements:
#' cell_types <- sort(unique(colData(sce_tr)[[column]]))
#' 
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom nnet multinom
#' @importFrom dplyr sample_n
train_nb <- function(x,y, cell_types) {
  
  flds <- caret::createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)
  x <- scale(x)
  
  metrics <- suppressWarnings({ 
    lapply(flds, function(test_idx) {
      # fit <- naive_bayes(x[-test_idx,], y[-test_idx])
      # p <- stats::predict(fit, newdata = x[test_idx,])
      df_train <- as.data.frame(x[-test_idx,])
      df_train$y <- y[-test_idx]
      df_train <- sample_n(df_train, min(5000, nrow(df_train)))
      
      df_test <- as.data.frame(x[test_idx,])
      df_test$y <- y[test_idx]
      
      fit <- multinom(y~., data = df_train)
      p <- stats::predict(fit, df_test)
      
      overall <- yardstick::bal_accuracy_vec(y[test_idx], p)
      scores <- sapply(cell_types, function(ct) {
        yardstick::f_meas_vec(factor(y[test_idx] == ct, levels=c("TRUE", "FALSE")), 
                factor(p == ct, levels = c("TRUE", "FALSE")))
      })
      tibble(
        what = c("Overall", as.character(cell_types)),
        score=c(overall, scores)
      )
  }) %>% bind_rows()
  })

  metrics
}

#' Make the expression and correlation heatmaps
#' 
#' @param sce A SingleCellExperiment object
#' @param markers A list storing three lists of markers: 
#' recommended_markers, scratch_markers, and top_markers
#' @param column A vector storing a column
#' @param display A vector storing selected heatmap display option
#' @param normalization A vector holding either 'Expression' or 'z-score'
#' @param pref_assay Assay loaded, default is 'logcounts'
#' 
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom stats cor
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly layout
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
  
  if(display == "Marker-marker correlation") {
    
    cc <- cor(t(mat))
    cor_map <- plot_ly(z=cc, 
            type='heatmap',
            x = rownames(cc),
            y = colnames(cc)) %>% 
      layout(title='Correlation')
    
    return(cor_map)
  } else {
    
    
    
    expression_map <- plot_ly(z=t(mat), 
            type='heatmap',
            x = rownames(mat),
            y = colnames(mat)) %>% 
      layout(title='Expression')

    return(expression_map)
  }
 
}

#' Suggest a set of redundant genes to remove
#' 
#' @param cmat A correlation matrix of *the full expression* calculated as
#' `expression <- as.matrix(assay(sce, pref_assay)[markers$top_markers,])`
#' `cmat <- cor(t(expression))`
#' @param n_genes Number of genes to suggest to remove
#' 
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
#' @param obj An R object storing metadata in a .rds file
#' 
#' TODO: accept anndata
#' 
read_input_scrnaseq <- function(sce) {

  if(methods::is(sce, 'Seurat')) {
    ## Convert to SingleCellExperiment
    sce <- Seurat::as.SingleCellExperiment(sce)
  } 
  stopifnot(methods::is(sce, 'SingleCellExperiment'))
  
  genes_to_remove <- grepl("^RP[L|S]|^MT-|^HSP|^FOS$|^JUN|MALAT1", rownames(sce))
  sce <- sce[!genes_to_remove,]
  sce
} 

#' Return the legend of the cell type colours
#' 
#' @param palette A vector storing hex colours
#' 
#' @importFrom ggplot2 theme_bw geom_histogram scale_fill_manual guides
get_legend <- function(palette) {
  ## Draw the ggplot
  df <- tibble(r = stats::runif(length(palette)),
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
#' 
#' @param gene_id An element in the marker list
#' @param df_antibody A dataframe storing antibody information
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
#' 
#' @param gene_id An element in the marker list
#' @param df A dataframe storing antibody information
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
#' 
#' @param gene_to_replace Gene to be replaced
#' @param sce A SingleCellExperiment object
#' @param pref_assay Assay loaded, default is 'logcounts'
#' @param n_correlations Number of markers to replace gene_to_replace
compute_alternatives <- function(gene_to_replace, sce, pref_assay = "logcounts", n_correlations) {
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

#' Parse out antibody applications
#' 
#' TODO: when the antibody application list is fixed, this can
#' be parsed ahead of time
get_antibody_applications <- function(antibody_info, 
                                      column_gene = 'Symbol',
                                      column_application = 'Listed Applications') {
  app_split <- strsplit(antibody_info[[ column_application ]], ',')
  app_split <- lapply(app_split, stringr::str_trim)
  names(app_split) <- antibody_info[[ column_gene ]]
  
  unique_applications <- sapply(app_split, `[`) %>% 
    unlist() %>% 
    table() %>% 
    sort(decreasing = TRUE)
  
  # Let's keep only those with > 500 genes/targets
  unique_applications <- unique_applications[unique_applications > 500]
  
  application_gene_map <- lapply(names(unique_applications), function(application) {
    v <- sapply(app_split, function(x) application %in% x)
    names(v[v])
  })
  names(application_gene_map) <- names(unique_applications)
  
  ## Convert to label format
  unique_applications2 <- names(unique_applications)
  names(unique_applications2) <- paste0(names(unique_applications), " (", unique_applications, ")")
  
  list(
    unique_applications = unique_applications2,
    application_gene_map = application_gene_map
  )
  
}

#' Picks out allowed genes corresponding
#' to antibody applications
get_allowed_genes <- function(selected_applications, applications_parsed, sce) {
  single_cell_genes <- rownames(sce)
  
  ## Get rid of RP[L|S] + MT + MALAT
  single_cell_genes <- single_cell_genes[!grepl("^RP[L|S]|^MT-|^MALAT", single_cell_genes)]
  
  antibody_genes <- unique(unlist(applications_parsed$application_gene_map))
  
  if(is.null(selected_applications)) {
    ## If no antibody application is selected,
    ## return all genes
    return(intersect(single_cell_genes, antibody_genes))
  }

  
  # ## Need to convert from the displayed application
  # ## name (which includes the number) back to the key
  # selected_applications <- plyr::mapvalues(selected_applications, 
  #                                          from = applications_parsed$unique_applications,
  #                                          to = names(applications_parsed$unique_applications))
  # 
  intersect(single_cell_genes, unique(unlist(applications_parsed$application_gene_map[selected_applications])))
}

#' Get the reactable for the modal dialogue when add
#' markers for a given cell type
#' @importFrom tibble rownames_to_column
get_cell_type_add_markers_reactable <- function(fm, current_markers) {
  fm <- fm[!rownames(fm) %in% current_markers,]
  fm <- fm[fm$summary.logFC > 0,]
  fm <- as.data.frame(head(fm, 50))
  fm <- rownames_to_column(fm, 'Gene')
  fm <- fm[,c("Gene", "FDR", "summary.logFC")]
  list(fm = fm,
       reactable = 
    reactable(
      fm,
      selection = "multiple",
      onClick = "select"
    )
)
}
