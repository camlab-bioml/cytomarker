#' Read the input single-cell RNA-seq dataset
#' from compressed RDS file. Accepts either:
#' (1) SingleCellExperiment
#' (2) Seurat object (converted to SingleCellExperiment)
#' 
#' @param sce_path Input uploaded path
#' @importFrom tools file_ext
#' @importFrom Matrix rowSums
read_input_scrnaseq <- function(sce_path) {
  
  library(SingleCellExperiment, quiet = T)
  library(SummarizedExperiment, quiet = T)
  library(zellkonverter, quiet = T)
  
  sce <- NULL ## object we're going to return
  
  if(file_ext(sce_path) == "h5ad") {
    ## We'll assume this is an h5ad file
    sce <- zellkonverter::readH5AD(sce_path)
  } else {
    ## We'll assume this is an rds
    sce <- readRDS(sce_path)
    if(!(isTruthy(methods::is(sce, 'SingleCellExperiment')) || isTruthy(methods::is(sce, 'Seurat')))) {
      sce <- NULL
    }
    
    if(methods::is(sce, 'Seurat')) {
      ## Convert to SingleCellExperiment
      sce <- Seurat::as.SingleCellExperiment(sce)
    } 
  }

  if (isTruthy(sce)) {
    if ("logcounts" %in% names(assays(sce))) {
      sce <- sce[, Matrix::colSums(logcounts(sce)) > 0]
    }
  }
  ## Remove cells with no reads - would cause issue with logNormCounts from scater
  sce
} 

#' This function detects whether the logcounts assay are indeed logcounts by finding
#' the residual of the logcount expression assays' rowsums divided by one. If this
#' is zero they must be counts and are converted to logcounts
#' @param sce Uploaded SingleCellExperiment
detect_assay_and_create_logcounts <- function(sce){
  library(scater, quiet = T)
  library(SummarizedExperiment, quiet = T)
  if ("logcounts" %in% names(assays(sce))) {
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
    sce <- scater::logNormCounts(sce)
  }
  
  sce <- sce[, Matrix::colSums(logcounts(sce)) > 0]
  
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
#' @param grch38 a human reference genome that is compatible with annotables
#' @param sce SingleCellExperiment object
#' @return a list with the first element being the status of the search, and the second being the genes found
check_rownames_for_ensembl<- function(sce, grch38){
  library(scuttle, quiet = T)
  ensemble_genes <- grepl("ENSG[0-9]*", rownames(sce))
  if(any(ensemble_genes)){
    sub_sce <- sce[ensemble_genes]
    rownames(sub_sce) <- gsub("\\.[0-9]", "", rownames(sub_sce))
    
    gene_map <- select(grch38, ensgene, symbol) |>
      filter(symbol != "") |>
      deframe()
    
    filtered_sce <- sub_sce[rownames(sub_sce) %in% names(gene_map)]
    
    rownames(filtered_sce) <- gene_map[rownames(filtered_sce)]
    
    filtered_mat <- scuttle::sumCountsAcrossFeatures(filtered_sce,
                                            rownames(filtered_sce),
                                            average = TRUE,
                                            assay.type = 'logcounts')
    
    filtered_sce <- SingleCellExperiment(list(logcounts = filtered_mat),
                                         colData = SingleCellExperiment::colData(sce))
    
    genes_found <- calculate_proportion_in_annotables(rownames(filtered_sce), grch38)
    # calculate the proportion as the number of retained ensembl vs. original gene identifiers
    genes_found$proportion <- nrow(filtered_sce) / nrow(sce)
    
    
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
#' @param grch38 Dataframe with two columns: \code{ensgene} and \code{symbol}, (originally from the annotables package)
#' @param remove_confounding_genes If true (default) genes that frequently confound single cell
#' analyses are removed. The following genes are removed: ribosobal proteins, mitochondrial ribosomal
#' and other mitochondrial proteins, heat shock proteins, Jun, Fos and Malat1
#'
#' @return genes A vector of gene names
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
    warning("Less than 50% of the genes in your dataset
                           match the database.")
  }else{
    ## STEP 2: Check if gene names are in rowData
    step2 <- check_rowData_for_hugo(sce, grch38)
    if(step2$status == "use_rowData_names"){
      clean_sce <- sce
      rownames(clean_sce) <- step2$genes
    }else if(step2$status == "low_number_of_gene_matches"){
      clean_sce <- sce
      rownames(clean_sce) <- step2$genes
      warning("Less than 50% of the genes in your dataset
                             match the database.")
    }else if(step1 != "rownames_are_null"){
      ## STEP 3: check if rownames are ensembl
      step3 <- check_rownames_for_ensembl(sce, grch38)
      
      if(step3$status == "use_ensembl_names"){
        clean_sce <- step3$sce
      }else if(step3$status == "low_number_of_gene_matches"){
        clean_sce <- step3$sce
        warning("Less than 50% of the genes in your dataset
                               match the database.")
      }else if(step3$status == "no_human_ensID"){
        stop("Detected ensembl ID's. However, these do not appear to be human.
                               Only human gene names are supported at the moment.
                               Please check the documentation for details.")
      }else if(step3$status == "did_not_find_ensembl"){
        stop("No human gene names were found in your dataset.
                               This could be because your genes are from a different species or
                               there are other mismatches. Please check the documentation for details.")
      }else{
        stop("No human gene names were found in your dataset.
                               This could be because your genes are from a different species or
                               there are other mismatches. Please check the documentation for details.")
      }
      
    }
  }
  if(exists("clean_sce")){
    clean_sce
  }else{
    stop("No human gene names were found in your dataset.
          This could be because your genes are from a different species or
          there are other mismatches. Please check the documentation for details.")
  }
}


#' Create a data frame table of counts for a heterogeneity category
#' @importFrom dplyr filter mutate rename arrange mutate_if
#' @param sce A SingleCellExperiment object
#' @param metadata_column the string name of a metadata column held within sce, on which to create frequency counts
create_table_of_hetero_cat <- function(sce, metadata_column) {
  table(sce[[metadata_column]]) |> 
    as.data.frame() |> 
    mutate(`Proportion Percentage` = 100*(Freq/sum(Freq))) |> 
    mutate_if(is.numeric, round, digits = 3) |>
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
  if (!is.null(pal_seed)) {
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
  sce$keep_for_analysis <- ifelse(sce[[input_column]] %in% vec_to_keep & sce$keep_for_analysis == "Yes",
                                  "Yes", "No")
  
  return(sce)
}

#' Given a vector of membership during sce sub-sampling, overwrite the keep_for_analysis column in the sce
#' @importFrom magrittr %>% 
#' @importFrom dplyr mutate
#' @param sce A SingleCellExperiment object
#' @param vec_to_keep A vector containing the values of an input column that should be retained for analysis
create_keep_vector_during_subsetting <- function(sce, vec_to_keep) {
  placeholder_frame <- data.frame(index = seq(1, nrow(SingleCellExperiment::colData(sce)))) %>%
    mutate(in_subsample = ifelse(index %in% vec_to_keep, "Yes", "No"))
  
  sce$in_subsample <- placeholder_frame$in_subsample
  sce$keep_for_analysis <- ifelse(sce$in_subsample == "Yes" & sce$keep_for_analysis == "Yes",
                                  "Yes", "No")
  sce$in_subsample <- NULL
  return(sce)
}


#' Update the keep_for_analysis vector in an SCE to exclude null and NA values
#' @param sce A SingleCellExperiment object
#' @param input_column The string name of the desired input column in sce
remove_null_and_va_from_cell_cat <- function(sce, input_column) {
  sce$keep_for_analysis <- ifelse((!is.null(sce[[input_column]]) |
                                     sce[[input_column]] != "NA" |
                                     sce[[input_column]] != "Na" |
                                     sce[[input_column]] != "na" |
                                     !is.na(sce[[input_column]])),
                                  "Yes", "No")
  
  return(sce)
}

#' Detect if the sce has any UMAP dimension assays 
detect_umap_dims_in_sce <- function(sce) {
  return(SingleCellExperiment::reducedDimNames(sce)[grepl("UMAP|umap|Umap|uMap|uMAP",
                          SingleCellExperiment::reducedDimNames(sce))])
}

#' Check if the majority of genes are not human in an SCE
#' @param sce The SingleCellExperiment object
#' @importFrom inferorg inferorg
check_for_human_genes <- function(sce) {
  inferorg_results <- inferorg(rownames(sce))
  if (inferorg_results$organism != "human") {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Remove confounding genes in the dataset 
#' @param gene_vector A vector of genes to clean
remove_confounding_genes <- function(gene_vector) {
  genes_to_remove <- grepl("^RP[L|S]|^MRP[L|S]|^MT-|^HSP|^FOS|^JUN|^MALAT1", gene_vector)
  clean_genes <- gene_vector[!genes_to_remove]
  return(clean_genes)
}
