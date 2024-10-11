
#' create a list of the run parameters for download or run logging
#' @importFrom yaml write_yaml
#' @importFrom zip zip
#' @importFrom htmlwidgets saveWidget
#' @param marker_list Vector of markers int the current panel
#' @param input_file Name of the scRNAseq file used to generate the panel
#' @param assay_used Name of scRNAseq assay used (common include counts or logcounts)
#' @param het_source Name of cell stratifying column in colData
#' @param panel_size Size of the panel created
#' @param cell_cutoff_value Minimum number of cells per cell type to include in analysis
#' @param subsample Whether subsetting of the dataset was done or not
#' @param subsample_number If subsetting was used, the number of cells used in analysis
#' @param marker_strat Type of panel building used. Options are cell-based or cell-free
#' @param antibody_apps List of antibody applications associated with the current panel
#' @param selected_cell_types All cell subtypes from het_source included in analysis
#' @param precomputed_umap_used If precompute UMAP coordinates were used or not
#' @param num_cells Number of cells in final analysis
#' @param num_genes Number of genes in final analysis
#' @param metrics Table of F1 scores per cell type from regression analysis
create_run_param_list <- function(
                          marker_list,
                          input_file,
                          assay_used,
                          het_source,
                          panel_size,
                          cell_cutoff_value,
                          subsample,
                          subsample_number,
                          marker_strat,
                          antibody_apps,
                          selected_cell_types,
                          precomputed_umap_used,
                          num_cells,
                          num_genes,
                          metrics) {
  
  ## Write a config list:
  config <- list(
    Time = as.character(Sys.time()),
    `Input file` = input_file,
    `Number of columns (cells)` = num_cells,
    `Number of rows (features)` = num_genes,
    `Assay used` = assay_used,
    `Heterogeneity source` = het_source,
    `Target panel size` = panel_size,
    `Min Cell Category cutoff` = cell_cutoff_value,
    `Subsampling Used` = subsample,
    `Subsampling number` = subsample_number,
    `Selected marker panel` = marker_list$top_markers,
    # `Scratch marker panel` = ifelse(is_empty(marker_list$scratch_markers),
    #                                 "None", marker_list$scratch_markers),
    `Scratch marker panel` = marker_list$scratch_markers,
    `Marker strategy` = marker_strat,
    # `Antibody applications` = ifelse(is_empty(antibody_apps),
    #                                  "None", antibody_apps),
    `Antibody applications` = antibody_apps,
    `Cell Types Analyzed` = selected_cell_types,
    `Pre-computed UMAP` = precomputed_umap_used,
    `Run Metrics` = metrics
  )
  
  return(config)
}


#' Download all marker data to a zip file
#' 
#' @importFrom yaml write_yaml
#' @importFrom zip zip
#' @importFrom htmlwidgets saveWidget
#' @importFrom readr write_tsv
#' @importFrom tibble tibble
#' @importFrom rmarkdown render
#' @param zip_filename Filepath where the zip output should be written
#' @param config Run configuration list in YAML format
#' @param plots List of session plots (UMAP, heatmap)
#' @param heatmap Heatmap plot to export
#' @param antibody_table Table of antibody products in the current panel
#' @param markdown_path Path to the package markdown report to export as HTML
#' @param run_metrics Table of F1 scores per cell type from regression analysis
#' @param overall_metrics Balanced accuracy score from F1 scores for the current run
#' @param markers_with_cell_type Named vector of current panel with associated cell type
download_data <- function(zip_filename,
                          config,
                          plots,
                          heatmap,
                          antibody_table,
                          markdown_path,
                          run_metrics,
                          overall_metrics,
                          markers_with_cell_type) {

    tmpdir <- tempdir()
    current_date <- Sys.Date()

    paths_zip <- list() # contains the filenames to save
    paths_report <- list() # contains the file names to be read in by report
    
    paths_zip$config <- file.path(tmpdir, paste0("config-", current_date, ".yml"))
    paths_report$config <- file.path(tmpdir, paste0("config-", current_date, ".tsv"))
    
    write_yaml(config, paths_zip$config)
  
    # selected_markers <- config$`Selected marker panel`
    
    config_df <- tibble::enframe(config) %>%
      dplyr::mutate(value = purrr::map_chr(value, toString)) |>
      `colnames<-`(c("Parameter", "Value")) |>
      filter(! Parameter %in% c("Input file", "Run Metrics"))
    
    param_list <- list(tmpdir = tmpdir,
                       df = antibody_table,
                       marker_selection = markers_with_cell_type,
                       heatmap = heatmap,
                       umap = list(all = plots$all, top = plots$top),
                       metric = list(plot = plots$metric,
                                     table = run_metrics),
                       overall = list(score = overall_metrics$score,
                                      counts = overall_metrics$counts),
                       config = config_df)
    
    paths_report$tmpdir <- paste0(tmpdir, "/")
    
    render(markdown_path, 
           output_file = paste0(paths_report$tmpdir, "report-", current_date, ".html"),
           output_dir = paths_report$tmpdir,
           params = param_list)
    
    paths_zip$report <- paste0(paths_report$tmpdir, "report-", current_date, ".html")
    zip(zipfile = zip_filename, files = unlist(paths_zip), mode = "cherry-pick") 
}


#' read back in a saved config yaml to resume analysis
#' @importFrom yaml read_yaml
#' @param yaml Path to a yaml file with a saved configuration
read_back_in_saved_yaml <- function(yaml) {
  return(read_yaml(yaml))
}

