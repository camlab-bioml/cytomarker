
#' create a list of the run parameters for download or run logging
#' @importFrom yaml write_yaml
#' @importFrom zip zip
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly subplot
create_run_param_list <- function(
                          marker_list,
                          input_file,
                          assay_used,
                          het_source,
                          panel_size,
                          cell_cutoff_value,
                          subsample,
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
#' @importFrom plotly subplot
download_data <- function(zip_filename,
                          config,
                          plots,
                          heatmap,
                          antibody_table) {

    tmpdir <- tempdir()
  
    paths_zip <- list()
    
    paths_zip$config <- file.path(tmpdir, paste0("config-", Sys.Date(), ".yml"))
    
    write_yaml(config, paths_zip$config)
    
    
    paths_zip$marker_selection <- file.path(tmpdir, paste0("markers-", Sys.Date(), ".txt"))
    
    selected_markers <- config$`Selected marker panel`
    write_lines(selected_markers, paths_zip$marker_selection)
    
    if (isTruthy(antibody_table)) {
      paths_zip$df <- file.path(tmpdir, paste0("Antibody-info-", Sys.Date(), ".tsv"))
      write.table(antibody_table, paths_zip$df, quote = F, row.names = F, sep = "\t")
    }    
    
    if (isTruthy(heatmap)) {
      paths_zip$heatmap <- file.path(tmpdir, paste0("Heatmap-", Sys.Date(), ".html"))
      saveWidget(heatmap, paths_zip$heatmap)
    }    
      
    if (isTruthy(plots$all_plot) & isTruthy(plots$top_plot)) {
        # paths_zip$umap <- file.path(tmpdir, paste0("UMAP-", Sys.Date(), ".pdf"))
        paths_zip$umap <- file.path(tmpdir, paste0("UMAP-", Sys.Date(), ".html"))
        umap_plt <- subplot(plots$all_plot, plots$top_plot) %>% 
          layout(title = 'Cytosel UMAP, all markers & top markers')
        saveWidget(umap_plt, paths_zip$umap)
        
      }
      
      if (isTruthy(plots$metric_plot)) {
        # paths_zip$metric <- file.path(tmpdir, paste0("metric-", Sys.Date(), ".pdf"))
        paths_zip$metric <- file.path(tmpdir, paste0("metrics-", Sys.Date(), ".html"))
        saveWidget(plots$metric_plot, paths_zip$metric)
      }
    
    zip(zipfile = zip_filename, files = unlist(paths_zip), mode = "cherry-pick") 
    if(file.exists(paste0(zip_filename, ".zip"))) {file.rename(paste0(zip_filename, ".zip"), zip_filename)}
}


#' read back in a saved config yaml to resume analysis
#' @importFrom yaml read_yaml
read_back_in_saved_yaml <- function(yaml) {
  return(read_yaml(yaml))
}

