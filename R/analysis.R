#' Compute the findMarkers outputs and store
#' 
#' Note: this used to be part of get_markers but needs
#' split off after colouring
#' 
#' @param sce A SingleCellExperiment object
#' @param columns A vector storing columns
#' @param pref_assay Assay loaded
#' @param allowed_genes List of genes to be included in the expression analysis
#' @param p_val_type The type of p-value to use in marker group calculation. Can be either all or any.
#' 
compute_fm <- function(sce, columns, pref_assay, allowed_genes,
                       p_val_type = "all") {
  
  library(scran)
  
  fms <- lapply(columns,
                  function(col) {
    
    test_type <- ifelse(pref_assay == "counts", "binom", "t")
    ag <- intersect(rownames(sce), allowed_genes)
    fm <- scran::findMarkers(sce[ag,], SummarizedExperiment::colData(sce)[[col]],
                      test.type = test_type, 
                      # BPPARAM = MulticoreParam(),
                      pval.type = p_val_type,
                      assay.type = pref_assay)
    
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
#' @param allowed_genes Set of allowed genes
#' @param in_session whether the function is being called in a shiny session or not
#' @importFrom dplyr mutate tally group_by filter pull slice_head arrange summarize ungroup
#' @import geneBasisR
get_markers <- function(fms, panel_size, marker_strategy, sce, allowed_genes,
                        in_session = T) {
  
  cell_types_wout_markers <- c()
  columns <- names(fms)
  original_multimarkers <- c()
  
  marker <- list(recommended_markers = c(), scratch_markers = c(), top_markers = c())
  
  if(marker_strategy == "geneBasis") {
    sce2 <- retain_informative_genes(sce[allowed_genes,], n = 10*panel_size)
    genes <- gene_search(sce2, n_genes=panel_size)
    marker <- list(
      recommended_markers = genes$gene[!is.na(genes$gene)],
      scratch_markers = c(),
      top_markers = genes$gene[!is.na(genes$gene)]
    )
  } else {
    for(col in columns) {
      fm <- fms[[col]]
      n <- length(fm)
      
      top_select <- ceiling((panel_size + 10) / n)
      
      recommended <- marker$recommended_markers
      recommended_df <- tibble(markers = c(), cell_type = c())
      top <- marker$top_markers
      scratch <- marker$scratch_markers
      
      highly_expressed <- tibble(markers = c(), cell_type = c())
      
      ## Create a vector of cell types for which no markers were found
      cell_types_wout_markers <- c()
      original_multimarkers <- c()
      for(i in seq_len(n)) {
        f <- fm[[i]]
        f <- f[!(rownames(f) %in% recommended),] |> as.data.frame()
        
        ## Only keep markers that are over-expressed
        f[is.na(f)] <- 0
        f <- f[f$summary.logFC > 0,]
        
        if(nrow(f) > 0){
          selected_markers <- rownames(f)[seq_len(top_select)]
          recommended_df <- bind_rows(recommended_df, 
                                      tibble(marker = selected_markers,
                                             cell_type = names(fm)[i],
                                             summary.logFC = f[selected_markers,]$summary.logFC))
          
          expressed_markers <- rownames(f)[seq_len(25)]
          highly_expressed <- bind_rows(highly_expressed, 
                                      tibble(marker = expressed_markers,
                                             cell_type = names(fm)[i],
                                             summary.logFC = f[expressed_markers,]$summary.logFC))
          
          
        }else{
          cell_types_wout_markers <- c(cell_types_wout_markers, names(fm)[i])
          # message(paste("No markers were found for the following cell types: ",
          #               names(fm)[i], sep = ""))
        }
        
        
        #recommended <- c(top, rownames(f)[seq_len(top_select)])
        #top <- c(top, rownames(f)[seq_len(top_select)])
      }
      
      initial_recommendations <- recommended_df
      
      # Iteratively remove markers until number of markers equals panel size
      # this needs to happen iteratively to ensure cell types have roughly
      # the same number of markers (otherwise if one cell type has all markers
      # with low lfc, all markers for that cell type could be removed)
      while(length(unique(recommended_df$marker)) > panel_size){
        # Find cell types that have the most number of markers. 
        # These are the ones markers can be removed from
        markers_per_cell_type <- group_by(recommended_df, cell_type) %>% 
          tally() %>% 
          mutate(max_markers = max(n)) %>% 
          filter(n == max_markers)
        
        # Select single marker to remove based on lowest lfc
        # sometimes the same marker is selected for several cell types
        # in this case the average lfc is taken
        remove <- recommended_df %>% 
          filter(cell_type %in% markers_per_cell_type$cell_type) %>% 
          group_by(marker) %>% 
          summarize(mean_lfc = mean(summary.logFC)) %>% 
          arrange(mean_lfc) %>% 
          slice_head(n = 1) %>% 
          pull(marker)
        
        recommended_df <- filter(recommended_df, marker != remove)
      }
      
      recommended <- unique(recommended_df$marker)
      
      # print cell types that have only markers that are expressed in multiple cell types
      genes_to_ignore <- c()
      count <- 0
      
      if (!is.null(recommended) & !is_empty(recommended) & length(allowed_genes) > 1000) {
        multimarkers <- create_multimarker_frame(initial_recommendations, recommended)
        
        genes_to_ignore <- c(genes_to_ignore, unique(multimarkers$marker))
        original_multimarkers <- c(original_multimarkers, unique(multimarkers$marker))
        
        while (nrow(multimarkers) > 0) {

          count <- count + 1
          genes_to_ignore <- c(genes_to_ignore, unique(multimarkers$marker))
          recommended <- recommended[!recommended %in% unique(multimarkers$marker) &
                                       !recommended %in% genes_to_ignore]
          multimarkers <- multimarkers |> group_by(cell_type) |> slice_head(n=1) |>
            ungroup() |> group_by(marker) |> slice_head(n=1)

          for (i in unique(multimarkers$cell_type)) {
            f <- fm[[i]]
            dont_use <- subset(highly_expressed, cell_type != i)
            f <- f[!rownames(f) %in% recommended &
                     !rownames(f) %in% multimarkers$marker &
                     !rownames(f) %in% genes_to_ignore &
                     !rownames(f) %in% dont_use$marker,] |> as.data.frame()

            ## Only keep markers that are over-expressed
            f[is.na(f)] <- 0
            f <- f[f$summary.logFC > 0,]
            new_markers_add <- rownames(f)[seq_len(subset(multimarkers, cell_type == i)$num_markers)]
            recommended <- c(recommended, new_markers_add)
          }
          multimarkers <- create_multimarker_frame(initial_recommendations, recommended)
        }
      }
      
      scratch <- unique(scratch)
      top <- recommended #unique(top)
    }
    marker <- list(recommended_markers = recommended[!is.na(recommended)],
                   scratch_markers = scratch[!is.na(scratch)],
                   top_markers = top[!is.na(top)])
  }
  return(list(marker = marker, missing = cell_types_wout_markers,
              multimarkers = original_multimarkers))
}


#' Collect possible multimarkers from the output of marker finding
#' @importFrom dplyr mutate tally group_by filter pull slice_head arrange summarize ungroup
#' @param frame A dataframe containing the markers selected in `findMarkers`
#' @param recommended_markers A vector of the currently recommended markers for a panel
create_multimarker_frame <- function(frame, recommended_markers) {
  return(frame |>
    filter(marker %in% recommended_markers) |> group_by(marker) |> 
    # get number of cell types each marker is expressed in
    mutate(cell_types_expressed_in = dplyr::n()) |> ungroup() |> 
    group_by(marker) |> slice_head(n=1) |> group_by(cell_type) |> 
    # get number of markers for a specific cell type
    mutate(num_markers = length(unique(marker)),
           highest = max(cell_types_expressed_in),
           lowest = min(cell_types_expressed_in)) |> ungroup() |>
    arrange(summary.logFC) |> filter(lowest > 1))
}


#' Given a list of top markers and the fms, get the associated
#' cell types
#' 
#' @param markers A list storing three lists of markers: 
#' recommended_markers, scratch_markers, and top_markers
#' @param fms Stored findMarkers outputs
get_associated_cell_types <- function(markers, fms) {
  fm <- fms[[1]] # For now, we're only doing this for the first
  
  all_markers <- unique(unlist(markers[c('recommended_markers', 
                        'scratch_markers', 'top_markers')]))
  
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
#' @param pref_assay Assay loaded
#' @param precomputed_vals Whether or not precomputed UMAP vals exist for the input SCE
#' @param dim_col If precomputed values exist, the name of the dimension holding the coordinates
#' @param only_top_markers Whether or not the UMAP computed will contain a subset of markers (Default is False)
#' @param markers_to_use The list of markers to compute the UMAP
#' 
#' @importFrom tibble tibble
get_umap <- function(sce, columns, pref_assay, precomputed_vals = NULL, dim_col = NULL,
                     only_top_markers = F, markers_to_use) {
  
  sce <- sce[rownames(sce) %in% markers_to_use,]
  
  # set max components to 25, or set components to the number of top markers if fewer than 25
  # num_comp_use <- ifelse(marker_num <= 25, marker_num, 25)
  
  if (!isTruthy(precomputed_vals) | isTRUE(only_top_markers)) {
    sce <- scater::runUMAP(sce, exprs_values = pref_assay,
                   # n_threads = availableCores(),
                   # ncomponents = num_comp_use, 
                   # pca = num_comp_use,
                   # BPPARAM = MulticoreParam(),
                   external_neighbors	= T)
    
    df <- as.data.frame(SingleCellExperiment::reducedDim(sce, 'UMAP')) |> `colnames<-`(c("UMAP_1", "UMAP_2"))
    
  } else {
    
    df <- as.data.frame(SingleCellExperiment::reducedDim(sce, dim_col)) |> `colnames<-`(c("UMAP_1", "UMAP_2"))
  }
    
  for(column in columns) {
    df[[column]] <- as.data.frame(SummarizedExperiment::colData(sce))[[column]]
    df[[column]] <- as.character(df[[column]])
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
#' @param pref_assay Assay loaded
get_scores <- function(sce, columns, mrkrs, pref_assay, max_cells = 5000) {
  
  max_cells <- min(ncol(sce), max_cells)
  sce_tr <- sce[mrkrs, sample(ncol(sce), max_cells, replace=FALSE)]
  
  scores <- list()
  
  for(column in columns) {
    scores[[column]] <- get_scores_one_column(sce_tr, column, mrkrs, pref_assay, max_cells)
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
#' @param pref_assay Assay loaded
get_scores_one_column <- function(sce_tr, column, mrkrs, pref_assay, max_cells = 5000) {
  
  x <- t(as.matrix(assay(sce_tr, pref_assay)))
  y <- factor(SummarizedExperiment::colData(sce_tr)[[column]])
  
  cell_types <- sort(unique(SummarizedExperiment::colData(sce_tr)[[column]]))
  
  train_nb(x,y, cell_types)
}

#' Train a ML model
#' 
#' Note this currently trains logistic regression
#' 
#' @param x A matrix calculated as `x <- t(as.matrix(assay(sce_tr)[[column]]))`
#' @param y A factor calculated as `y <- factor(SummarizedExperiment::colData(sce_tr)[[column]])`
#' @param cell_types A sorted vector with unique elements:
#' cell_types <- sort(unique(SummarizedExperiment::colData(sce_tr)[[column]]))
#' 
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom nnet multinom
#' @importFrom dplyr sample_n
train_nb <- function(x,y, cell_types) {

    flds <- caret::createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)
    x <- scale(x)
    x <- x[,colMeans(is.na(x)) == 0]
    
    metrics <- lapply(flds, function(test_idx) {
        # fit <- naive_bayes(x[-test_idx,], y[-test_idx])
        # p <- stats::predict(fit, newdata = x[test_idx,])
        df_train <- as.data.frame(x[-test_idx,])
        df_train$y <- y[-test_idx]
        df_train <- sample_n(df_train, min(5000, nrow(df_train)))
        
        # write.csv(df_train, "df_train.csv", quote = F, na = "NA")
        
        df_test <- as.data.frame(x[test_idx,])
        df_test$y <- y[test_idx]
        
        fit <- tryCatch({
          multinom(y~., data = df_train, trace = FALSE, MaxNWts = 100000)
        },
        error = function(cond) {
          return(NULL)
        }, 
        warning = function(cond) {
          return(NULL)
        }
        )
        
        if (isTruthy(fit)) {
          p <- stats::predict(fit, df_test)
          
          overall <- yardstick::bal_accuracy_vec(y[test_idx], p)
          scores <- sapply(cell_types, function(ct) {
            yardstick::f_meas_vec(factor(as.character(y[test_idx]) == as.character(ct), levels=c("TRUE", "FALSE")), 
                                  factor(as.character(p) == as.character(ct), levels = c("TRUE", "FALSE")))
          })
          tibble(what = c("Overall", as.character(cell_types)),
            score=c(overall, scores)
          )
    } else {
      tibble(what = c("Overall", as.character(cell_types)),
        score=c(0, rep(0, length(as.character(cell_types)))))
    }
        
    }) %>% bind_rows()
        # fit <- glmnet::glmnet(
        #   x = iris.x, y = iris.y,
        #   family = "multinomial",
        #   lambda = 0,
        #   type.multinomial = "grouped"
        # )
    
  return(metrics)
  
  
}

# train_2 <- function(x, y, cell_types) {
#   y <- df_train$y
#   x <- dplyr::select(df_train, -y)
#   cvfit <- cv.glmnet(as.matrix(x), y, family = "multinomial")
#   
# }


#' Suggest a set of redundant genes to remove
#' 
#' @param cmat A correlation matrix of *the full expression* calculated as
#' `expression <- as.matrix(assay(sce, pref_assay)[markers$top_markers,])`
#' `cmat <- cor(t(expression))`
#' @param n_genes Number of genes to suggest to remove
#' 
suggest_genes_to_remove <- function(cmat, n_genes=10) {
  rg <- c()
  
  for(i in seq_len(n_genes)) {
    lgl <- !(rownames(cmat) %in% rg)
    fcs <- caret::findCorrelation(cmat[lgl, lgl], cutoff = 0.01)
    gene_to_remove <- colnames(cmat[lgl, lgl])[ fcs[1] ]
    rg <- c(rg, gene_to_remove)
  }
  
  rg
  
  
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
#' @param pref_assay Assay loaded
#' @param n_correlations Number of markers to replace gene_to_replace
#' @param allowed_genes A list of genes permitted to be suggested. Prevents suggestion of unwanted genes such as mitochondrial, etc.
#' @param fms The output of findMarkers
#' @importFrom dplyr arrange
compute_alternatives <- function(gene_to_replace, sce, pref_assay, n_correlations,
                                 allowed_genes, fms) {
  x <- as.matrix(assay(sce, pref_assay))[,sample(ncol(sce), min(5000, ncol(sce)))]
  
  y <- x[gene_to_replace, ]
  
  yo <- x[rownames(x) != gene_to_replace,]
  
  correlations <- cor(t(yo), y)
  
  alternatives <- data.frame(Gene = rownames(yo), Correlation = round(correlations[,1], 3))
  alternatives <- alternatives[alternatives$Gene %in% allowed_genes,]
  alternatives <- alternatives[!is.na(alternatives$Correlation),]
  alternatives <- alternatives[order(-(alternatives$Correlation)),]
  alternatives <- alternatives[1:n_correlations,]
  associated <- get_associated_cell_types(list(recommended_markers = alternatives$Gene), fms)
  alternatives <- merge(alternatives, data.frame(Gene = names(associated),
                                                `Association` = associated),
                        by = "Gene") |> arrange(-Correlation)
  rownames(alternatives) <- NULL
  
  alternatives
  
}

#' Get the reactable for the modal dialogue when add
#' markers for a given cell type
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate_if mutate_at
#' @param fm Output from findMarkers
#' @param current_markers vector of markers in the current panel
get_cell_type_add_markers_reactable <- function(fm, current_markers) {
  fm <- fm[!rownames(fm) %in% current_markers,]
  fm <- fm[fm$summary.logFC > 0,]
  fm <- as.data.frame(head(fm, 50))
  fm <- rownames_to_column(fm, 'Gene')
  fm <- fm[,c("Gene", "FDR", "summary.logFC")]
  fm$FDR <- signif(fm$FDR, 3)
  fm$summary.logFC <- round(fm$summary.logFC, digits = 3)
  list(fm = fm,
       reactable = 
         reactable(
           fm,
           selection = "multiple",
           onClick = "select"
         )
  )
}

