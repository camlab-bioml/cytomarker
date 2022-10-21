
#' Download all marker data to a zip file
#' 
#' @importFrom yaml write_yaml
#' @importFrom zip zip
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly subplot
#' @importFrom readr write_tsv
#' @importFrom tibble tibble
download_data <- function(zip_filename,
                          current_markers,
                          plots,
                          heatmap,
                          input_file,
                          assay_used,
                          het_source,
                          panel_size,
                          cell_cutoff_value,
                          subsample,
                          antibody_table,
                          marker_strat,
                          antibody_apps,
                          selected_cell_types,
                          precomputed_umap_used,
                          num_cells,
                          num_genes) {
    browser()
    tmpdir <- tempdir()
  
    paths_zip <- list()
    paths_quarto <- list()
    
    paths_zip$config <- file.path(tmpdir, paste0("config-", Sys.Date(), ".yml"))
    paths_quarto$config <- file.path(tmpdir, paste0("config-", Sys.Date(), ".tsv"))
    
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
      `Selected marker panel` = current_markers$top_markers,
      `Scratch marker panel` = current_markers$scratch_markers,
      `Marker strategy` = marker_strat,
      `Antibody applications` = antibody_apps,
      `User selected cells` = selected_cell_types,
      `Pre-computed UMAP` = precomputed_umap_used
    )
    
    config_df <- tibble(
      p = names(config),
      v = lapply(config, function(x) if(length(x) > 1) paste(x, collapse = ", ") else x) |> unlist(use.names = FALSE)
    )

    write_yaml(config, paths_zip$config)
    write_tsv(config_df, paths_quarto$config)
    
    
    paths_quarto$marker_selection <- file.path(tmpdir, paste0("markers-", Sys.Date(), ".txt"))
    
    selected_markers <- current_markers$top_markers
    write_lines(selected_markers, paths_quarto$marker_selection)
    
    if (isTruthy(antibody_table)) {
      paths_quarto$df <- file.path(tmpdir, paste0("Antibody-info-", Sys.Date(), ".tsv"))
      write.table(antibody_table, paths_quarto$df, quote = F, row.names = F, sep = "\t")
    }    
    
    if (isTruthy(heatmap)) {
      paths_quarto$heatmap <- file.path(tmpdir, paste0("Heatmap-", Sys.Date(), ".rds"))
      saveRDS(heatmap, paths_quarto$heatmap)
    }    
      
    if (isTruthy(plots$all_plot) & isTruthy(plots$top_plot)) {
        # paths_zip$umap <- file.path(tmpdir, paste0("UMAP-", Sys.Date(), ".pdf"))
        paths_quarto$umap <- file.path(tmpdir, paste0("UMAP-", Sys.Date(), ".rds"))
        umap_plt <- subplot(plots$all_plot, plots$top_plot) %>% 
          layout(title = 'Cytosel UMAP, all markers & top markers')
        saveRDS(umap_plt, paths_quarto$umap)
      }
      
      if (isTruthy(plots$metric_plot)) {
        # paths_zip$metric <- file.path(tmpdir, paste0("metric-", Sys.Date(), ".pdf"))
        paths_quarto$metric <- file.path(tmpdir, paste0("metrics-", Sys.Date(), ".rds"))
        saveRDS(plots$metric_plot, paths_quarto$metric)
      }
    
    quarto::quarto_render("report/testing-quarto.qmd", output_format = "html",
                          execute_params = paths_quarto)
    zip(zipfile = zip_filename, files = unlist(paths_zip), mode = "cherry-pick") 
    if(file.exists(paste0(zip_filename, ".zip"))) {file.rename(paste0(zip_filename, ".zip"), zip_filename)}
}


#' read back in a saved config yaml to resume analysis
#' @importFrom yaml read_yaml
read_back_in_saved_yaml <- function(yaml) {
  return(read_yaml(yaml))
}

