
#' create a list of the run parameters for download or run logging
#' @importFrom yaml write_yaml
#' @importFrom zip zip
#' @importFrom htmlwidgets saveWidget
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
download_data <- function(zip_filename,
                          config,
                          plots,
                          heatmap,
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
read_back_in_saved_yaml <- function(yaml) {
  return(read_yaml(yaml))
}

