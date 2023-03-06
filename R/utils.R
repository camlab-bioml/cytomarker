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


round3 <- function(x) format(round(x, 1), nsmall = 3)

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

#' Convert a SingleCellExperiment column to a factor or a character
#' @param sce a SingleCellExperiment object
#' @param input_column a metadata column in the sce input
convert_column_to_character_or_factor <- function(sce, input_column) {
  if (is.numeric(sce[[input_column]])) {
    sce[[input_column]] <- as.character(factor(sce[[input_column]],
                                               levels = sort(unique(sce[[input_column]]))
    ))
  } else {
    sce[[input_column]] <- as.character(sce[[input_column]])
  }
  sce
}

#' #' replace a null, NA, or empty value with a string representation "None
#' #' @param val the list element to convert
#' replace_na_null_empty <- function(val) {
#'   if (is_empty(val)) {
#'     return ("None")
#'   } else if (is.na(val)) {
#'     return ("None")
#'   } else {
#'     return(val)
#'   }
#' }

#' Filter out antibody table rows by a vector list of antibody applications
#' @param applications A vector of antibody applications
#' @param row The string in a row representing the applications for a specific gene in the table
in_row <- function(applications, row) {
  return (length(intersect(applications, unlist(strsplit(row, ", ")))) > 0)
}


#' Get a list of aliases for the marker list (all possible gene IDs that match to an Entrez ID)
#' @param marker_list A vector of markers in gene symbol format
#' @param gene_mapper A data frame that maps gene symbols to Entrez IDs
#' @param merge Boolean whether the aliases should be returned along with the original markers
get_gene_aliases <- function(marker_list, gene_mapper, allowed_genes) {
  
  # get the list of entrez ids for the uploaded panel
  possible_entrez <- gene_mapper[gene_mapper$alias_symbol %in% marker_list,]
  # expand the possible symbols to include symbols with the same entrez id
  all_possible_symbols <- gene_mapper[gene_mapper$gene_id %in% possible_entrez$gene_id,]
  
  # filter the new
  all_possible_symbols <- all_possible_symbols[all_possible_symbols$alias_symbol %in% allowed_genes,]
  
  # keep only genes from the original upload that are allowed
  filtered_upload <- marker_list[marker_list %in% allowed_genes]
  
  # get a list of genes that were replaced with an alias
  replaced <- unique(all_possible_symbols$alias_symbol[!all_possible_symbols$alias_symbol %in% marker_list])
  
  reduced_table <- merge(all_possible_symbols |> `colnames<-`(c("entrez_id", "Alias")),
                          possible_entrez |> `colnames<-`(c("entrez_id", "Original Name")),
                          by = "entrez_id") |> filter(Alias != `Original Name`)
  
  # get a list of the original marker list that was replaced 
  # do not include symbols that could not be mapped to Exntra Ids as they will be kept
  removed <- marker_list[!marker_list %in% all_possible_symbols$alias_symbol &
                           marker_list %in% possible_entrez$alias_symbol]
  
  list(merged = unique(c(filtered_upload, all_possible_symbols$alias_symbol)),
                            original = marker_list, removed = removed, replaced = replaced,
       table = reduced_table)
  
}





