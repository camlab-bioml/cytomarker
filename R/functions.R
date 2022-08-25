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
#' @param pref_assay Assay loaded
#' @importFrom SingleCellExperiment colData
#' @importFrom scran findMarkers
#' @importFrom parallel mclapply
#' @importFrom parallelly availableCores
compute_fm <- function(sce, columns, pref_assay, allowed_genes) {

  fms <- mclapply(columns, mc.cores = availableCores(), function(col) {
    
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
#' @param allowed_genes Set of allowed genes
#' @import geneBasisR
#' @importFrom dplyr mutate tally group_by filter pull slice_head arrange summarize
get_markers <- function(fms, panel_size, marker_strategy, sce, allowed_genes) {
  
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
        
      top_select <- ceiling((panel_size + 10) / n)
      
      recommended <- marker$recommended_markers
      recommended_df <- tibble(markers = c(), cell_type = c())
      top <- marker$top_markers
      scratch <- marker$scratch_markers
        
      ## Create a vector of cell types for which no markers were found
      cell_types_wout_markers <- c()
      for(i in seq_len(n)) {
        f <- fm[[i]]
        f <- f[!(rownames(f) %in% recommended),]
        
        ## Only keep markers that are over-expressed
        f[is.na(f)] <- 0
        f <- f[f$summary.logFC > 0,]
        
        if(nrow(f) > 0){
          selected_markers <- rownames(f)[seq_len(top_select)]
          recommended_df <- bind_rows(recommended_df, 
                                      tibble(marker = selected_markers,
                                             cell_type = names(fm)[i],
                                             summary.logFC = f[selected_markers,]$summary.logFC))
        }else{
          cell_types_wout_markers <- c(cell_types_wout_markers, names(fm)[i])
        }
        #recommended <- c(top, rownames(f)[seq_len(top_select)])
        #top <- c(top, rownames(f)[seq_len(top_select)])
      }
      
      if(!is.null(cell_types_wout_markers)){
        throw_error_or_warning(type = 'notification',
                               message = paste("No markers were found for the following cell types: ",
                                               paste(cell_types_wout_markers, 
                                                     collapse = ", "), 
                                               ". This is likely because there are too few cells of these types."),
                               duration = 10)
      }
      
      #recommended_df <- group_by(recommended_df, marker, cell_type)
      
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
      scratch <- unique(scratch)
      top <- recommended #unique(top)
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
set_current_markers_safely <- function(markers, fms, default_type = NULL) {
  
  if (is.null(markers$associated_cell_types)) {
    markers$associated_cell_types <- get_associated_cell_types(markers, fms)
  }
  if (is.list(markers$associated_cell_types)) {
    markers$associated_cell_types <- unlist(markers$associated_cell_types)
  }
    
  markers
}

#' Compute the UMAP 
#'
#' @param sce A SingleCellExperiment object
#' @param columns A vector storing columns
#' @param pref_assay Assay loaded
#' 
#' @importFrom scater runUMAP
#' @importFrom tibble tibble
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom parallelly availableCores
#' @importFrom BiocParallel MulticoreParam
get_umap <- function(sce, columns, pref_assay) {
  
  sce <- runUMAP(sce, exprs_values = pref_assay,
                 n_threads = availableCores(),
                 ncomponents = 20, 
                 pca = 20,
                 # BPPARAM = MulticoreParam(workers = availableCores()),
                 external_neighbors	= T)
  
  df <- as.data.frame(tibble(
    UMAP1 = as.numeric(reducedDim(sce, 'UMAP')[,1]),
    UMAP2 = as.numeric(reducedDim(sce, 'UMAP')[,2])
  ))
  
  for(column in columns) {
    df[[column]] <- as.data.frame(colData(sce))[[column]]
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
  y <- factor(colData(sce_tr)[[column]])
  
  cell_types <- sort(unique(colData(sce_tr)[[column]]))
  
  train_nb(x,y, cell_types)
}

#' Train a ML model
#' 
#' Note this currently trains logistic regression
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
#' @importFrom parallel mclapply
#' @importFrom parallelly availableCores
train_nb <- function(x,y, cell_types) {
  
  flds <- caret::createFolds(y, k = 10, list = TRUE, returnTrain = FALSE)
  x <- scale(x)
  
  metrics <- suppressWarnings({ 
    mclapply(flds, mc.cores = availableCores(), function(test_idx) {
      # fit <- naive_bayes(x[-test_idx,], y[-test_idx])
      # p <- stats::predict(fit, newdata = x[test_idx,])
      df_train <- as.data.frame(x[-test_idx,])
      df_train$y <- y[-test_idx]
      df_train <- sample_n(df_train, min(5000, nrow(df_train)))
      
      df_test <- as.data.frame(x[test_idx,])
      df_test$y <- y[test_idx]
      
      fit <- multinom(y~., data = df_train, trace = FALSE)
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

# train_2 <- function(x, y, cell_types) {
#   y <- df_train$y
#   x <- dplyr::select(df_train, -y)
#   cvfit <- cv.glmnet(as.matrix(x), y, family = "multinomial")
#   
# }

#' Make the expression and correlation heatmaps
#' 
#' @param sce A SingleCellExperiment object
#' @param markers A list storing three lists of markers: 
#' recommended_markers, scratch_markers, and top_markers
#' @param column A vector storing a column
#' @param display A vector storing selected heatmap display option
#' @param normalization A vector holding either 'Expression' or 'z-score'
#' @param pref_assay Assay loaded
#' 
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom stats cor
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly layout
create_heatmap <- function(sce, markers, column, display, normalization, pref_assay) {

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
#' @param sce_path Input uploaded path
#' @importFrom tools file_ext
#' @importFrom zellkonverter readH5AD
#' @importFrom SingleCellExperiment logcounts
#' @importFrom Matrix rowSums
read_input_scrnaseq <- function(sce_path) {
  sce <- NULL ## object we're going to return
  
  if(file_ext(sce_path) == "h5ad") {
    ## We'll assume this is an h5ad file
    sce <- readH5AD(sce_path)
  } else {
    ## We'll assume this is an rds
    sce <- readRDS(sce_path)
    if(!(isTruthy(methods::is(sce, 'SingleCellExperiment')) || isTruthy(methods::is(sce, 'Seurat')))) {
      invalid_modal()
    }
      
    if(methods::is(sce, 'Seurat')) {
      ## Convert to SingleCellExperiment
      sce <- Seurat::as.SingleCellExperiment(sce)
    } 
  }
  ## Remove cells with no reads - would cause issue with logNormCounts from scater
  sce <- sce[, Matrix::colSums(logcounts(sce)) > 0]
  sce
} 

#' This function detects whether the logcounts assay are indeed logcounts by finding
#' the residual of the logcount expression assays' rowsums divided by one. If this
#' is zero they must be counts and are converted to logcounts
#' @param sce Uploaded SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom scater logNormCounts
detect_assay_and_create_logcounts <- function(sce){
  count_sums <- assays(sce)$logcounts |> 
    rowSums()
  ## calculate residuals after dividing by 1
  modulo_residuals <- lapply(count_sums, function(x) x %% 1) |> 
    unlist() |> 
    sum()
  
  ## If divisible by one must be counts, this converts to logcounts
  if(modulo_residuals == 0){
    assays(sce)[['counts']] <- assays(sce)$logcounts
    assays(sce)[['logcounts']] <- NULL
    sce <- logNormCounts(sce)
  }
  
  sce
}

##### [ PARSE GENE NAME FUNCTIONS ] #####

#' Helper function for `parse_gene_names`
#' Calculates the overlap of gene vector with annotables
#' 
#' @importFrom magrittr %>% 
#' @importFrom dplyr mutate filter pull
calculate_proportion_in_annotables <- function(g, annotation){
  # Count how many of the input genes were found in annotables
  genes_found <- factor(g %in% annotation$symbol, levels = c(TRUE, FALSE))
  genes_found <- table(genes_found) %>% 
    as.data.frame() %>% 
    mutate(total = sum(Freq),
           prop_found = Freq / total)
  
  prop_found <- filter(genes_found, genes_found == TRUE) %>% 
    pull(prop_found)
  genes_found <- filter(genes_found, genes_found == TRUE) %>% 
    pull(Freq)
  
  list(proportion = prop_found, gene_num = genes_found)
}

#' Checks if sce rownames can be used
check_rownames_for_hugo <- function(sce, grch38){
  if(is.null(rownames(sce))){
    error_status <- "rownames_are_null"
  }else{
    rowname_genes <- rownames(sce)
    
    # Calculate overlap with annotables
    genes_found <- calculate_proportion_in_annotables(rowname_genes, grch38)
    prop_found <- genes_found$proportion
    genes_found <- genes_found$gene_num
    
    # Infer organism
    inferred_type <- inferorg(rowname_genes)
    
    # Save as gene_names if conditions match
    if(prop_found > 0.5 && inferred_type$organism == 'human'){
      error_status <- "use_rownames"
    }else if(prop_found <= 0.5 && genes_found > 100 && inferred_type$organism == 'human'){
      error_status <- "low_number_of_gene_matches"
    }else{
      if(inferred_type$organism == 'human'){
        error_status <- "rownames_do_not_match_criteria"
      }else{
        error_status <- "did_not_find_hugo_in_rownames"
      }
    }
  }
  
  error_status
}

#' Checks if there are any columns in rowData that contain gene names
#' @importFrom stats na.omit
check_rowData_for_hugo <- function(sce, grch38){
  # List of possible rowData column names
  possible_symbol_rowData_names <- c("Gene", "gene", "Genes", "genes", "geneID", "GeneID", 
                                     "symbol", "Symbol", "Gene_ID", "gene_id", "gene_ID",
                                     "gene_symbol", "gene_Symbol", "geneSymbol", "GeneSymbol",
                                     "genesymbol", "genessymbol", "genessymbols", "genesymbols")
  
  # This finds all the ID's of rowData(sce) that match possible_symbol_rowData_names
  gene_match_idx <- na.omit(match(possible_symbol_rowData_names, colnames(rowData(sce))))
  #gene_match_idx <- na.omit(match(colnames(rowData(sce)), possible_symbol_rowData_names))
  
  if(length(gene_match_idx) > 0){
    # Go through every match and find the one with the highest score
    colData_matches <- lapply(gene_match_idx, function(x){
      g <- rowData(sce)[[x]]
      inferred_type <- inferorg(g)
      
      tibble(organism = inferred_type$organism,
             format = inferred_type$format,
             confidence_format = inferred_type$confidence_format,
             colData_column_idx = x)
    }) %>% bind_rows() %>% 
      filter(organism == "human" & format == "symbol") %>% 
      arrange(-confidence_format)
    
    # Check if there are any rows left after filtering down to human and symbol
    if(nrow(colData_matches) > 0){
      # Take top hit by format confidence score
      best_symbol_idx <- colData_matches$colData_column_idx[1]
      colData_genes <- rowData(sce)[[best_symbol_idx]]
      
      genes_found <- calculate_proportion_in_annotables(colData_genes, grch38)
      if(genes_found$proportion > 0.5){
        output <- list(status = "use_rowData_names",
                       genes = colData_genes)
      }else if(genes_found$proportion <= 0.5 && genes_found$gene_num > 100){
        output <- list(status = "low_number_of_gene_matches",
                       genes = colData_genes)
      }else{
        output <- list(status = "rowdata_does_not_match_criteria", genes = NULL)
      }
    }else{
      output <- list(status = "no_symbols_found_in_rowdata", genes = NULL)
    }
  }else{
    output <- list(status = "no_symbols_found_in_rowdata", genes = NULL)
  }
  
  output
}

#' Checks if the rownames of the sce contain human ensembl ID's
#' @importFrom dplyr select
#' @importFrom scuttle sumCountsAcrossFeatures
#' @importFrom SingleCellExperiment SingleCellExperiment
check_rownames_for_ensembl<- function(sce, grch38){
  ensemble_genes <- grepl("ENSG[0-9]*", rownames(sce))
  if(any(ensemble_genes)){
    sce <- sce[ensemble_genes]
    rownames(sce) <- gsub("\\.[0-9]", "", rownames(sce))
    
    gene_map <- select(grch38, ensgene, symbol) %>% 
      filter(symbol != "") %>% 
      deframe()
    
    filtered_sce <- sce[rownames(sce) %in% names(gene_map)]
    
    rownames(filtered_sce) <- gene_map[rownames(filtered_sce)]
    
    filtered_mat <- sumCountsAcrossFeatures(filtered_sce,
                                            rownames(filtered_sce),
                                            average = TRUE,
                                            assay.type = 'logcounts')
    
    filtered_sce <- SingleCellExperiment(list(logcounts = filtered_mat))
    colData(filtered_sce) <- colData(sce)
    
    genes_found <- calculate_proportion_in_annotables(rownames(filtered_sce), grch38)
    
    if(genes_found$proportion > 0.5){
      output <- list(status = "use_ensembl_names",
                     sce = filtered_sce)
    }else if(genes_found$proportion <= 0.5 && genes_found$gene_num > 100){
      output <- list(status = "low_number_of_gene_matches",
                     sce = filtered_sce)
    }else{
      output <- list(status = "ensembl_does_not_match_criteria", sce = NULL)
    }
    
    output
  }else if(any(grepl("ENS", rownames(sce)))){
    list(status = "no_human_ensID", sce = NULL)
  }else{
    list(status = "did_not_find_ensembl", sce = NULL)
  }
}

#' Parses the gene names in a sce object
#' This function tries to find gene symbols by going through several steps:
#' 1. Count the number and proportion of rownames that are listed in the
#'    symbol column of the annotables grch38 object. 
#'    If the proportion is more than 50% and the symbols are human they are used
#'    If the proportion is less than 50%, more than 100 genes are found and 
#'    the symbols are human they are still used
#' 2. Check if there are any human gene names in `rowData(sce)` and if more
#'    than 50%/100 genes are in the symbol column of the annotables grch38 object
#' 3. If rownames are not NULL, check if they are ensembl ID's and convert to symbol
#' 
#' @param sce SingleCellExperiment object
#' 
#' @return genes A vector of gene names, or should it return an sce?
#' 
#' @importFrom inferorg inferorg
#' @importFrom dplyr mutate filter pull bind_rows
#' @importFrom tibble tibble deframe
#' @importFrom magrittr %>% 
parse_gene_names <- function(sce, grch38){
  ## STEP 1: Check if rownames can be used
  step1 <- check_rownames_for_hugo(sce, grch38)
  
  if(step1 == "use_rownames"){
    clean_sce <- sce
  }else if(step1 == "low_number_of_gene_matches"){
    clean_sce <- sce
    throw_error_or_warning(type = 'notification',
                           message = "Less than 50% of the genes in your dataset
                           match the database.",
                           notificationType = 'warning')
  }else{
    ## STEP 2: Check if gene names are in rowData
    step2 <- check_rowData_for_hugo(sce, grch38)
    if(step2$status == "use_rowData_names"){
      clean_sce <- sce
      rownames(clean_sce) <- step2$genes
    }else if(step2$status == "low_number_of_gene_matches"){
      clean_sce <- sce
      rownames(clean_sce) <- step2$genes
      throw_error_or_warning(type = 'notification',
                             message = "Less than 50% of the genes in your dataset
                             match the database.",
                             notificationType = 'warning')
    }else if(step1 != "rownames_are_null"){
      ## STEP 3: check if rownames are ensembl
      step3 <- check_rownames_for_ensembl(sce, grch38)
      
      if(step3$status == "use_ensembl_names"){
        clean_sce <- step3$sce
      }else if(step3$status == "low_number_of_gene_matches"){
        clean_sce <- step3$sce
        throw_error_or_warning(type = 'notification',
                               message = "Less than 50% of the genes in your dataset
                               match the database.",
                               notificationType = 'warning')
      }else if(step3$status == "no_human_ensID"){
        throw_error_or_warning(type = 'error',
                               message = "Detected ensembl ID's. However, these do not appear to be human.
                               Only human gene names are supported at the moment. 
                               Please check the documentation for details.")
      }else if(step3$status == "ensembl_does_not_match_criteria"){
        throw_error_or_warning(type = 'error',
                               message = "No human gene names were found in your dataset.
                               This could be because your genes are from a different species or
                               there are other mismatches. Please check the documentation for details.")
      }else{
        throw_error_or_warning(type = 'error',
                               message = "No human gene names were found in your dataset.
                               This could be because your genes are from a different species or
                               there are other mismatches. Please check the documentation for details.")
      }
      
    }
  if(exists("clean_sce")){
    # Remove RPL/S Mitochondrial, HSP, FOS, JUN and MALAT1 genes
    genes_to_remove <- grepl("^RP[L|S]|^MT-|^HSP|^FOS$|^JUN|MALAT1", rownames(clean_sce))
    clean_sce[!genes_to_remove,]
  }else{
    throw_error_or_warning(type = 'error',
                           message = "No human gene names were found in your dataset.
                           This could be because your genes are from a different species or
                           there are other mismatches. Please check the documentation for details.")
    }
  }
}

#' Global function to throw errors/warnings
#' 
#' @param type needs to be `error` or `notification`. `error` will create a popup error
#'  while `notification` will display a notification in the bottom right. Defaults to `notification`.
#' @param message the message to be displayed to the user
#' @param duration number of seconds to display message before it disappears. 
#' Defaults to NULL (does not close until closed by user).
#' @param notificationType Controls the colour of the message. Allowed values are 'default' (grey), 
#' 'message' (blue), 'warning' (yellow) and 'error' (red). Only applicable to notifications
#' @param errorTitle Title to display on the error message. Defaults to 'Error'
throw_error_or_warning <- function(type = 'notification', 
                                   message, 
                                   duration = NULL,
                                   notificationType = 'default',
                                   errorTitle = "Error"){
  
  # Change duration format so it can also be accepted by shinyalert
  # shinyalert requires timer = 0 if notification should be closed by user
  # and duration is in miliseconds as opposed to seconds for showNotification
  if(type == 'error'){
    if(is.null(duration)){
      duration <- 0
    }else{
      duration <- duration * 1000
    }
  }
  
  if(type == 'error'){
    shinyalert(
      title = errorTitle,
      text = message,
      type = "error",
      timer = duration,
      showConfirmButton = TRUE,
      confirmButtonCol = "#337AB7"
    )
  }else if(type == 'notification'){
    showNotification(
      message,
      type = notificationType,
      duration = duration)
  }
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
#' whether the gene is initially suggested or selected by the user
#' 
#' @param gene_id An element in the marker list
#' @param markers the list of markers
map_gene_name_to_antibody_icon <- function(gene_id, markers) {
  
  if(gene_id %in% markers$recommended_markers){
    return(icon("star"))
  }else{
    return(NULL)
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
#' @param pref_assay Assay loaded
#' @param n_correlations Number of markers to replace gene_to_replace
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


#' Create a data frame table of counts for a heterogeneity category
#' @importFrom dplyr filter mutate rename arrange
#' @param sce A SingleCellExperiment object
#' @param metadata_column the string name of a metadata column held within sce, on which to create frequency counts
create_table_of_hetero_cat <- function(sce, metadata_column) {
  table(sce[[metadata_column]]) |> 
    as.data.frame() |> 
    mutate(`Proportion Percentage` = 100*(Freq/sum(Freq))) |> 
    rename(!!metadata_column := "Var1") |>
    arrange(-Freq)
}


#' Create a vector of metadata values that pass a minimum count
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @param grouped_frame a data frame generated using `create_table_of_hetero_cat`
#' @param sce A SingleCellExperiment object
#' @param metadata_column The string name of a metadata column held within sce. Should be the same as the column used in grouped_frame.
#' @param min_counts the minimum counts to retain a category in grouped frame. Suggested minimum is 2. 
remove_cell_types_by_min_counts <- function(grouped_frame, sce, metadata_column, min_counts = 2) {
  
  keep_frame <- grouped_frame %>% filter(Freq >= min_counts)
  
  # keep_frame[,] <- lapply(df, function(x) {as.numeric(as.character(x))})
  
  return(as.vector(keep_frame[metadata_column][,1]))
}

#' Create the global cytosel palette with or without seeding. The palette begins with 12 color-blind friendly colours
#' then moves into 74 uniquely generated colors form brewer.pal, ending with 2 repeats from the first vector for
#' a final vector of 100 unique colours.. 
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @param pal_seed random seed used to shuffle the palette (Default is NULL)
create_global_colour_palette <- function(pal_seed = NULL) {
  
  # color blind friendly palette: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
  
  safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  unique_palette <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  if (! is.null(pal_seed)) {
    set.seed(pal_seed)
    unique_palette <- sample(unique_palette, length(unique_palette))
  }
  unique_palette <- unique_palette[!unique_palette %in% safe_colorblind_palette]
  
  return(c(safe_colorblind_palette, unique_palette, safe_colorblind_palette,
           unique_palette[1:2]))
}


#' Convert a text colour into black or white based on the RGB values of its background to improve visibility
#' @importFrom grDevices col2rgb
#' @param text_colour the colour to be converted into black or white
set_text_colour_based_on_background <- function(text_colour) {
  rgb_code <- as.vector(col2rgb(text_colour))
  # https://stackoverflow.com/questions/3942878/how-to-decide-font-color-in-white-or-black-depending-on-background-color
  if (rgb_code[1]*0.299 + rgb_code[2]*0.587 + rgb_code[3]*0.114 > 150) {
    return("#000000")
  } else {
    return("#ffffff")
  }
}

#' Given a vector of input values, create an sce column for retaining during analysis based on membership of the input column in the retention vector
#' @param sce A SingleCellExperiment object
#' @param vec_to_keep A vector containing the values of an input column that should be retained for analysis
#' @param input_column The string name of the desired input column in sce. The values in this column should have overlap with vec_to_keep
create_sce_column_for_analysis <- function(sce, vec_to_keep, input_column) {
  sce$keep_for_analysis <- ifelse(sce[[input_column]] %in% vec_to_keep,
                                  "Yes", "No")
  
  return(sce)
}

#' Given a vector of membership during sce sub-sampling, overwrite the keep_for_analysis column in the sce
#' @importFrom magrittr %>% 
#' @importFrom dplyr mutate
#' @param sce A SingleCellExperiment object
#' @param vec_to_keep A vector containing the values of an input column that should be retained for analysis
create_keep_vector_during_subsetting <- function(sce, vec_to_keep) {
  placeholder_frame <- data.frame(index = seq(1, nrow(colData(sce)))) %>%
    mutate(in_subsample = ifelse(index %in% vec_to_keep, "Yes", "No"))
  
  sce$in_subsample <- placeholder_frame$in_subsample
  sce$keep_for_analysis <- ifelse(sce$in_subsample == "Yes" & sce$keep_for_analysis == "Yes",
                                 "Yes", "No")
  sce$in_subsample <- NULL
  return(sce)
}
