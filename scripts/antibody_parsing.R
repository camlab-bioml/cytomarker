# These are deprecated function used for vendor-specific antibody catalogs
# to collect the antibody applications associated with a product

#' Parse out antibody applications
#'
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