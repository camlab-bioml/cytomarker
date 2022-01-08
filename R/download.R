
#' Download all marker data to a zip file
#' 
#' @importFrom yaml write_yaml
#' 
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
    
    paths_zip$marker_selection <-file.path(tmpdir, paste0("markers-", Sys.Date(), ".txt"))
    paths_zip$umap <- file.path(tmpdir, paste0("UMAP-", Sys.Date(), ".pdf"))
    paths_zip$heatmap <- file.path(tmpdir, paste0("heatmap-", Sys.Date(), ".pdf"))
    paths_zip$metric <- file.path(tmpdir, paste0("metric-", Sys.Date(), ".pdf"))
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
    
    selected_markers <- current_markers$top_markers
    write_lines(selected_markers, paths_zip$marker_selection)
    
    umap_plt <- cowplot::plot_grid(plots$all_plot, plots$top_plot, nrow=1)
    ggsave(paths_zip$umap, umap_plt, width=10, height=4)
    
    pdf(paths_zip$heatmap, onefile = TRUE)
    draw(heatmap)
    dev.off()
    
    ggsave(paths_zip$metric, plots$metric_plot, width=10, height=4)
    
    zip(zipfile = zip_filename, files = unlist(paths_zip), mode = "cherry-pick") 
    if(file.exists(paste0(zip_filename, ".zip"))) {file.rename(paste0(zip_filename, ".zip"), zip_filename)}
}

