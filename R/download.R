
#' Download all marker data to a zip file
#' 
#' @importFrom yaml write_yaml
#' @importFrom zip zip
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly subplot
download_data <- function(zip_filename,
                          current_markers,
                          plots,
                          heatmap,
                          input_file,
                          assay_used,
                          het_source,
                          panel_size) {

    tmpdir <- tempdir()
  
    paths_zip <- list()
    
    paths_zip$config <- file.path(tmpdir, paste0("config-", Sys.Date(), ".yml"))
    
    ## Write a config list:
    config <- list(
      Time = as.character(Sys.time()),
      `Input file` = input_file,
      `Assay used` = assay_used,
      `Heterogeneity source` = het_source,
      `Target panel size` = panel_size
    )
    write_yaml(config, paths_zip$config)
    
    
    paths_zip$marker_selection <-file.path(tmpdir, paste0("markers-", Sys.Date(), ".txt"))
    
    selected_markers <- current_markers$top_markers
    write_lines(selected_markers, paths_zip$marker_selection)
    
    if (isTruthy(heatmap)) {
      paths_zip$heatmap <-file.path(tmpdir, paste0("Heatmap-", Sys.Date(), ".html"))
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

