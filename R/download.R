
#' Download all marker data to a zip file
#' 
#' @importFrom yaml write_yaml
#' @importFrom zip zip
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly subplot
#' @importFrom readr write_tsv
#' @importFrom tibble tibble
#' @importFrom rmarkdown render
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
    tmpdir <- tempdir()
    current_date <- Sys.Date()
  
    paths_zip <- list() # contains the filenames to save
    paths_report <- list() # contains the file names to be read in by report
    
    paths_zip$config <- file.path(tmpdir, paste0("config-", current_date, ".yml"))
    paths_report$config <- file.path(tmpdir, paste0("config-", current_date, ".tsv"))
    
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
    
    config <- lapply(config, FUN = function(X) replace_na_null_empty(X))
    
    config_df <- tibble(
      Parameter = names(config),
      Value = lapply(config, function(x) if(length(x) > 1) paste(x, collapse = ", ") else x) |> 
        unlist(use.names = FALSE)
    )
    
    write_yaml(config, paths_zip$config)
    write_tsv(config_df, paths_report$config)
    
    paths_report$marker_selection <- file.path(tmpdir, paste0("markers-", current_date, ".txt"))
    
    selected_markers <- current_markers$top_markers
    write_lines(selected_markers, paths_report$marker_selection)
    
    if (isTruthy(antibody_table)) {
      paths_report$df <- file.path(tmpdir, paste0("Antibody-info-", current_date, ".tsv"))
      write.table(antibody_table, paths_report$df, quote = F, row.names = F, sep = "\t")
    }    
    
    if (isTruthy(heatmap)) {
      paths_report$heatmap <- file.path(tmpdir, paste0("Heatmap-", current_date, ".rds"))
      saveRDS(heatmap, paths_report$heatmap)
    }    
      
    if (isTruthy(plots$all_plot) & isTruthy(plots$top_plot)) {
        paths_report$umap <- file.path(tmpdir, paste0("UMAP-", current_date, ".rds"))
        umap_plt <- subplot(plots$all_plot, plots$top_plot) %>%
          layout(title = 'Cytosel UMAP, all markers & top markers')
        saveRDS(umap_plt, paths_report$umap)
      }
      
      if (isTruthy(plots$metric_plot)) {
        paths_report$metric <- file.path(tmpdir, paste0("metrics-", current_date, ".rds"))
        saveRDS(plots$metric_plot, paths_report$metric)
      }
    
    paths_report$tmpdir <- paste0(tmpdir, "/")
    render("report/rmarkdown-report.Rmd", 
           output_file = paste0(paths_report$tmpdir, "report-", current_date, ".html"),
           output_dir = paths_report$tmpdir,
           params = paths_report)
    
    paths_zip$report <- paste0(paths_report$tmpdir, "report-", current_date, ".html")
    zip(zipfile = zip_filename, files = unlist(paths_zip), mode = "cherry-pick") 
    if(file.exists(paste0(zip_filename, ".zip"))) {file.rename(paste0(zip_filename, ".zip"), zip_filename)}
}


#' read back in a saved config yaml to resume analysis
#' @importFrom yaml read_yaml
read_back_in_saved_yaml <- function(yaml) {
  return(read_yaml(yaml))
}

