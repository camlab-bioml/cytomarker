
context("Test basic Shiny app server functionality")

test_that("Server has basic functionality", {
  
  testServer(cytosel::cytosel(), expr = {
    
    # Verify the input sce and default + selected assays
    expect_equal(pref_assay(), "logcounts")
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")))
    
    expect_is(sce(), 'SingleCellExperiment')
    expect_equivalent(dim(sce()), c(2167, 100))
    session$setInputs(assay_select = "counts", assay = "counts")
    expect_equal(pref_assay(), "counts")
    
    session$setInputs(precomputed_dim = T, select_precomputed_umap = "UMAP",
                      possible_precomputed_dims = reducedDimNames(sce()))
    expect_true(use_precomputed_umap())
    expect_equal(umap_precomputed_col(), "UMAP")
    
    session$setInputs(coldata_column = "seurat_annotations",
                      min_category_count = 2,
                      show_cat_table = T)
    
    expect_equal(output$cell_cat_preview,
              "9 groupings in selected category, including B, Undetermined, CD14+ Mono and others")
    
    # expect all 9 cell types to be in the tally and selected category
    expect_equal(nrow(summary_cat_tally()), 9)
    expect_equal(length(specific_cell_types_selected()), 9)
    expect_equal(specific_cell_types_selected(), unique(sce()[["seurat_annotations"]]))                
    
    # Set analysis parameters that will not proceed
    session$setInputs(user_selected_cells = NULL,
                      panel_size = 32, min_category_count = 0,
                      subset_number = 99,
                      subsample_sce = T,
                      start_analysis = T,
                      add_selected_to_analysis = T,
                      marker_strategy = "fm")
    
    # if the min count is set below 2, do not proceed with analysis
    expect_false(proceed_with_analysis())
    
    # do not have a plot generated if the analysis cannot proceed
    expect_null(plots$top_plot)
    
    session$setInputs(min_category_count = 50,
                      start_analysis = T)
    
    expect_false(any_cells_present())
    
    expect_null(previous_metrics())
    
    # Change to passable input parameters
    session$setInputs(user_selected_cells = c("CD8 T", "Memory CD4 T",
                                              "Naive CD4 T",
                                              "Platelet"),
                      add_selected_to_analysis = T)
    
    
    session$setInputs(min_category_count = 2,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression",
                      tabs = NULL,
                      metrics_toggle = NULL,
                      select_aa = NULL,
                      panel_sorter = "Group by cell type",
                      start_analysis = T)
    
    # resetting the min count to 2 allows to proceed with analysis
    expect_true(proceed_with_analysis())
    expect_true(any_cells_present())
    
    # platelet does not meet the cutoff, so of the 9 cell types, 8 pass high enough
    expect_equal(length(cell_types_high_enough()), 8)
    expect_equal(c("CD8 T", "Memory CD4 T",
                   "Naive CD4 T"), cell_types_to_keep())
    
    # Ensure that only yes and no are in the keep_for_analysis column
    expect_true("Yes" %in% unique(sce()$keep_for_analysis))
    expect_true("No" %in% unique(sce()$keep_for_analysis))
    expect_equal(length(unique(sce()$keep_for_analysis)), 2)
    
    # Ensure that the dimensions od the sce are not reduced after subsetting
    expect_equivalent(dim(sce()), c(2167, 100))
    
    # do not keep Platelet even though it was not selected for removal
    expect_equivalent(subset(colData(sce()), 
                             seurat_annotations == "Platelet")$keep_for_analysis,
                      "No")
    
    # expect that all cell types have markers for them
    
    expect_null(cell_types_missing_markers())
    
    # Verify proper format of the findmarkers output and marker lists
    expect_equal("seurat_annotations", column())
    expect_equal(length(names(fms()[[1]])), 3)
    
    expect_equal(length(current_markers()), 3)
    
    expect_false(is.null(current_markers()$recommended_markers))
    expect_false(is.null(current_markers()$top_markers))
    
    # expect that the palette is the same length as the number of cells retained
    expect_equal(names(cytosel_palette()), c("CD8 T", "Memory CD4 T", "Naive CD4 T"))

    # prove that the marker-marker correlation plot was made
    
    expect_equal(class(heatmap())[1], "plotly")
    expect_equal(class(heatmap())[2], "htmlwidget")
    
    expect_equal(heatmap()[["x"]][["attrs"]][[1]][["x"]],
                 heatmap()[["x"]][["attrs"]][[1]][["y"]])
    
    # expect of 200 markers to get somewhere between 1 to 200 for redundancy
    expect_gt(length(as.character(heatmap()[["x"]][["attrs"]][[1]][["x"]])), 1)
    expect_lte(length(as.character(heatmap()[["x"]][["attrs"]][[1]][["x"]])),
                 200)
    
    # expect that the marker columns are set properly and the lengths are verified
    
    # expect that the marker columns are set properly and the lengths are verified
    session$setInputs(bl_top = current_markers()$top_markers[1:10],
                      bl_scratch = current_markers()$top_markers[11:15],
                      start_analysis = T)
    
    expect_equal(num_markers_in_selected(), 10)
    expect_equal(num_markers_in_scratch(), 5)
    
    expect_equal(output$scratch_marker_counts, "<B> Scratch Markers: 5 </B>")
    expect_equal(output$selected_marker_counts, "<B> Selected Markers: 10 </B>")
    
    expect_null(marker_sort())
    
    session$setInputs(panel_sorter = "Sort alphabetically")
    expect_equal(marker_sort(), "Sort alphabetically")
    
    session$setInputs(panel_sorter = "Group by cell type")
    expect_equal(marker_sort(), "Group by cell type")
    
    session$setInputs(tabs = "Metrics", metrics_toggle = "Current Run Metrics")
    
    expect_true(!is.null(current_overall_score()))
    
    # suggest markers to remove
    
    session$setInputs(suggest_gene_removal = T,
                      n_genes = 3)
  
    expect_equal(length(suggestions()), 3)
    
    # move 5 of the top markers into scratch and assert new numbers
    session$setInputs(markers_to_remove = suggestions(),
                      remove_suggested = T)
    
    # a gene that doesn't exist will not produce suggested replacements
    expect_equal(length(current_markers()$top_markers), 7)
    
    # generate suggested alternative markers
    session$setInputs(input_gene = "FAKE_GENE",
                      number_correlations = 10,
                      enter_gene = T)
    
    expect_null(replacements())
    
    # generate suggested alternative markers
    session$setInputs(input_gene = current_markers()$top_markers[3],
                      number_correlations = 10,
                      enter_gene = T)
    
    # replacements should be of length 10
    expect_false(is.null(output$alternative_markers))
    
    expect_equal(nrow(replacements()), 10)
    expect_equal(ncol(replacements()), 3)
    
    expect_true("CD8 T" %in% unique(replacements()$Association))
    
    expect_equal(0, nrow(replacements()[grepl("^RP[L|S]|^MT-|^MALAT",
                                              replacements()$Gene),]))
    
    session$setInputs(alternative_markers_rows_selected = c(1, 3, 4, 7, 8, 10),
                      send_markers = T,
                      send = T)
    
    expect_equal(length(current_markers()$top_markers), 13)
    
    # if you try to add a fake marker, the length will stay the same
    session$setInputs(markers_change_modal = T,
                      add_markers = "FAKE_GENE",
                      enter_marker = T)
    
    expect_equal(num_markers_in_selected(), 13)
    
    
    # if you try to add a marker that's already there, the length will stay the same
    session$setInputs(markers_change_modal = T,
                      add_markers = current_markers()$top_markers[1],
                      enter_marker = T)
    
    expect_equal(num_markers_in_selected(), 13)
    
    # add a gene that isn't in either selected or scratch to increase the count by 1
    different <- allowed_genes()[!allowed_genes() %in% current_markers()$top_markers]
    different <- different[!different %in% current_markers()$scratch_markers]
    
    session$setInputs(markers_change_modal = T,
                      add_markers = different[1],
                      enter_marker = T)
    
    expect_equal(num_markers_in_selected(), 14)
    
    # expect null since no markers have been suggested
    expect_true(is.null(selected_cell_type_markers()))
    
    # get 50 suggested markers for the specific cell type
    session$setInputs(markers_change_modal = T,
                      cell_type_markers = "CD8 T",
                      add_cell_type_markers = T)
    
    # expect the number of suggestions to be 50
    expect_equal(nrow(marker_suggestions()), 50)
    
    # expect no rows when looking for unwanted genes
    expect_equal(0, nrow(marker_suggestions()[grepl("^RP[L|S]|^MT-|^MALAT",
                                            marker_suggestions()$Gene),]))
    
    # expect the suggestion output table is created
    expect_false(is.null(output$cell_type_marker_reactable))
    session$setInputs(add_select_markers = T)
    
    # Simulate uploading own markers
    session$setInputs(uploadMarkers = list(datapath =
                                             test_path("upload_markers.txt")),
                     add_to_selected = T)
    
    # after adding 5 markers to 14, assert how many top markers there should be
    expect_equal(length(current_markers()$top_markers), 19)
    
    # if you try to upload a gene that doesn't exist, the length will remain the same
    session$setInputs(uploadMarkers = list(datapath =
                                             test_path("fake_upload.txt")),
                      add_to_selected = T)
    
    # after adding 5 markers to 18, assert how many top markers there should be
    expect_equal(length(current_markers()$top_markers), 19)
    
    # if you try to replace with a gene that doesn't exist, the length will remain the same
    session$setInputs(uploadMarkers = list(datapath =
                                             test_path("fake_upload.txt")),
                      replace_selected = T)
    
    # after adding 5 markers to 18, assert how many top markers there should be
    expect_equal(length(current_markers()$top_markers), 19)
    
    # Look for the uploaded markers in the top markers
    expect_true("EEF2" %in% current_markers()$top_markers)
    expect_true("TRAM1" %in% current_markers()$top_markers)
    
    # instead of add, replace markers
    session$setInputs(uploadMarkers = list(datapath =
                                             test_path("upload_markers.txt")),
                      replace_selected = T)
    
    # after replacing all top markers with 5, assert how many top markers there should be
    expect_equal(length(current_markers()$top_markers), 5)
    
    expect_false("LTB" %in% current_markers()$top_markers)
    
    copy_markers <- current_markers()
    
    expect_false(is.null(current_markers()))
    expect_equal(copy_markers, current_markers())
    
    session$setInputs(markers_to_remove = NULL,
                      remove_suggested= T)
    
    expect_equal(length(current_markers()$top_markers), 5)
    
    
    withr::with_tempdir({
      
      # skip_on_ci()
      session$setInputs(downloadData = T)
      # expect_true(file.exists(file.path(tempdir(), paste0("config-", Sys.Date(), ".yml"))))
      expect_false(is.null(output$downloadData))
    })
    
  })
})

test_that("Pre-setting the input rank lists persists in the current markers", {
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      assay = "counts", coldata_column = "seurat_annotations")
    
    session$setInputs(subsample_sce = T,
                      panel_size = 20,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression",
                      marker_strategy = "fm",
                      select_aa = NULL,
                      bl_top = c("EEF2", "RBM3", "TRAM1", "MSN", "FTL"),
                      bl_recommended = c("EEF2", "RBM3", "TRAM1", "MSN", "FTL"),
                      bl_scratch = c("FGR", "ACAP1"),
                      start_analysis = T)
    
    expect_equal(length(current_markers()$top_markers), length(input$bl_top))
    expect_equal(length(current_markers()$scratch_markers), length(input$bl_scratch))
    
  })
  
})


context("Test that Shiny app server can detect single assay")

test_that("Server can detect sce with only one assay", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_equal(pref_assay(), "logcounts")
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_one_assay.rds")),
                      min_category_count = 2,
                      panel_size = 32)
    
    # detect if only one assay is uploaded
    expect_equal(length(names(assays(sce()))), 1)
    
  })
  
})

gc()

context("Test re-upload and reset Shiny app server functionality")

#### Re-upload analysis ######
test_that("Re-upload works on server", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_false(reupload_analysis())
  
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      read_back_analysis = list(datapath =
                                                  test_path("test_config.yml")))
    # verify re-upload analysis reactive worked
    expect_true(reupload_analysis())
    
    # verify that the reactive values were populated from the yml
    expect_equal(length(specific_cell_types_selected()), 5)
    expect_equal(cell_min_threshold(), 10)
    expect_equal(length(markers_reupload()$top_markers), 18)
    expect_equal(length(markers_reupload()$scratch_markers), 6)
    expect_equal(length(specific_cell_types_selected()), 5)
    
    session$setInputs(panel_size = 200, coldata_column = "seurat_annotations",
                      subsample_sce = F,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression",
                      marker_strategy = "fm")
    
    session$setInputs(tabs = NULL,
                      metrics_toggle = NULL,
                      select_aa = NULL,
                      panel_sorter = "Group by cell type")
    
    session$setInputs(start_analysis = T)
    
    expect_equal(length(current_markers()$top_markers),
                 length(markers_reupload()$top_markers))
    
    expect_equal(length(current_markers()$scratch_markers),
                 length(markers_reupload()$scratch_markers))
    
  })
})

gc()

context("Test that resetting the cytosel environment works as expected")


#### Reset analysis ######
test_that("Reset works on server", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_false(reupload_analysis())
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      read_back_analysis = list(datapath =
                                                  test_path("test_config.yml")))
    # verify reupload analysis reactive worked
    expect_true(reupload_analysis())
    
    # verify that the reactive values were populated from the yml
    expect_equal(length(specific_cell_types_selected()), 5)
    expect_equal(cell_min_threshold(), 10)
    expect_equal(length(markers_reupload()$top_markers), 18)
    expect_equal(length(markers_reupload()$scratch_markers), 6)
    
    session$setInputs(bl_scratch = markers_reupload()$scratch_markers,
                      bl_top = markers_reupload()$top_markers)
    
    session$setInputs(create_reset = T, reset_marker_panel = T)
    session$setInputs(dismiss_marker_reset = T)
    
    # resetting the panel sets current markers and input rank lists to 0
    expect_true(reset_panel())
    expect_null(current_markers()$top_markers)
    expect_equal(num_markers_in_selected(), 0)
    expect_equal(num_markers_in_scratch(), 0)
    
    expect_true(proceed_with_analysis())
    
  })
})

context("Current panel with different genes throws error")

test_that("Error from current panel with different genes", {
  testServer(cytosel::cytosel(), expr = {
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      panel_size = 20,
                      coldata_column = "seurat_annotations",
                      min_category_count = 2)
    
    expect_true(valid_existing_panel())
    
    session$setInputs(bl_top = c("GENE_1", "GENE_2"),
                      start_analysis = T)
    
    # check that genes that are not in the rownames of the sce invalidate the panel
    expect_false(valid_existing_panel())
    
  })
})

gc()

context("Test UMAP, violin, and heatmap colouring changes")

test_that("Changing the UMAP, violin, and heatmap colourings work", {
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      assay = "counts", coldata_column = "seurat_annotations")
    
    session$setInputs(show_cat_table = T)
    expect_equal(cell_min_threshold(), 2)
    expect_equal(length(specific_cell_types_selected()),
                 length(unique(sce()[[input$coldata_column]])))
    
    session$setInputs(subsample_sce = F,
                      panel_size = 200,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression")
    
    session$setInputs(panel_sorter = "Group by cell type",
                      marker_strategy = "fm",
                      tabs = NULL,
                      metrics_toggle = NULL,
                      select_aa = NULL)
    
    session$setInputs(start_analysis = T)
    
    expect_equal(cell_min_threshold(), 2)
    expect_equal(length(specific_cell_types_selected()),
                 length(unique(sce()[[input$coldata_column]])))
    
    # expect only the first run outputs to be rendered 
    
    expect_false(is.null(output$current_run_name))
    expect_false(is.null(output$summary_run_current))
    expect_false(is.null(output$metrics_run_current))

    expect_null(previous_run_log())
    
    heatmap_1 <- heatmap()
    
    expect_true(is.null(umap_all()))
    expect_true(is.null(umap_top()))
  
    session$setInputs(tabs = "UMAP", umap_options = "Cell Type",
                      umap_panel_options = "S100A9", umap_panel_cols = T, 
                      show_umap_legend = T)
    
    expect_false(is.null(umap_all()))
    expect_false(is.null(umap_top()))
    
    # check defaults for UMAP plots
    expect_false(umap_top_gene())
    expect_false(umap_all_gene())
    expect_equal(umap_colouring(), "Cell Type")
    
    session$setInputs(tabs = "UMAP", umap_options = "Panel Marker",
    umap_panel_options = "S100A9", umap_panel_cols = T, show_umap_legend = T)
    
    # switching the UMAP to gene works
    expect_false(isFALSE(umap_top_gene()))
    expect_false(isFALSE(umap_all_gene()))
    expect_equal(umap_colouring(), "Panel Marker")
    
    viol_markers <- c("EEF2", "RBM3", "TRAM1", "MSN", "FTL")
    
    # setting the violin plots with genes works
    session$setInputs(genes_for_violin = viol_markers,
                      add_violin_genes = T, viol_viewer = "By Marker")
    expect_false(is.null(output$expression_violin))
    expect_true(all(viol_markers %in% current_markers()$top_markers))
    
    current_length <- length(current_markers()$top_markers)
    
    viol_1 <- output$expression_violin
    
    session$setInputs(genes_for_violin = viol_markers,
                      add_violin_genes = T, viol_viewer = "By Cell Type")
    
    expect_equal(length(current_markers()$top_markers), current_length)
    
    expect_false(identical(viol_1, output$expression_violin))
    
    # if set to NULL, still renders an empty violin plot
    session$setInputs(genes_for_violin = NULL)
    expect_false(is.null(output$expression_violin))
    
    # switching the heatmap category makes a different heatmap
    session$setInputs(display_options = "seurat_annotations",
    heatmap_expression_norm = "Expression", start_analysis = T)
    
    expect_false(identical(heatmap_1, heatmap()))
    
    session$setInputs(start_analysis = T)

    expect_true(isTruthy(previous_run_log()))
    expect_false(is.null(output$previous_run_1))
    expect_false(is.null(output$summary_prev_1))
    expect_false(is.null(output$metrics_run_prev_1))
    
  })
})

# if (file.exists(file.path(tempdir(), "/seurat_pbmc.rds"))) {
#   command <- paste('rm ', tempdir(), "/seurat_pbmc.rds", sep = "")
#   system(command)
# }

gc()

context("Test the loading of the curated datasets from dropbox")

test_that("Picking the curated dataset works as intended", {
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs(curated_dataset = T, curated_options = "Kidney",
                      coldata_column = "cell_ontology_class",
                      curated_compartments = c("Endothelial", "Immune", "Stromal"),
                      pick_curated = T)
    
    expect_equal(curated_selection(), "Kidney")
    
    expect_true(file.exists(file.path(tempdir(), "/Kidney.rds")))
    expect_equivalent(dim(sce()), c(23863, 750))
    
    expect_false("epithelial" %in% unique(sce()$compartment))
    expect_false(is.null(output$curated_set_preview))
    
    session$setInputs(select_precomputed_umap = "UMAP",
                      possible_precomputed_dims = reducedDimNames(sce()))
    
    # session$setInputs(start_analysis = T)
    
    # if reupload the same dataset, invalidate until verify if resetting or not
    
    session$setInputs(curated_dataset = T, curated_options = "Liver",
                      subset_number = 2000,
                      coldata_column = "cell_ontology_class",
                      curated_compartments = NULL,
                      pick_curated = T)
    
    expect_equal(curated_selection(), "Liver")
    
  })
  
})

test_that("Setting null compartments retains the full dataset", {
  testServer(cytosel::cytosel(), expr = {

    session$setInputs(curated_dataset = T, curated_options = "Kidney",
                  subset_number = 2000,
                  coldata_column = "cell_ontology_class",
                  curated_compartments = NULL,
                  pick_curated = T)
    expect_true(file.exists(file.path(tempdir(), "/Kidney.rds")))
    expect_equivalent(dim(sce()), c(23863, 881))

    expect_true("epithelial" %in% unique(sce()$compartment))
    
    expect_true(proceed_with_analysis())
    
    session$setInputs(panel_size = 0, start_analysis = T)
    
    expect_false(proceed_with_analysis())
    

})
  
})

test_that("Having an existing panel will warn for a reset on upload", {
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs( bl_top = c("EEF2", "RBM3", "TRAM1", "MSN", "FTL"),
                       bl_recommended = c("EEF2", "RBM3", "TRAM1", "MSN", "FTL"),
                       bl_scratch = c("GNLY", "FTL"))
    
    session$setInputs(curated_dataset = T, curated_options = "Kidney",
                      subset_number = 2000,
                      coldata_column = "cell_ontology_class",
                      curated_compartments = NULL,
                      pick_curated = T)
    expect_true(file.exists(file.path(tempdir(), "/Kidney.rds")))
    expect_equivalent(dim(sce()), c(23863, 881))
    
    expect_false(proceed_with_analysis())
    expect_true("epithelial" %in% unique(sce()$compartment))
    
    session$setInputs(coldata_column = "n_genes", assay = "logcounts",
                      precomputed_dim = "UMAP",
                      panel_size = 24)
    
    session$setInputs(start_analysis = T)
    
    expect_false(proceed_with_analysis())
    
  })
  
})

context("Test that finding markers with very few genes produces an error")

test_that("datasets with few genes produce errors on marker finding", {
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_few_genes.rds")),
                      assay = "logcounts", coldata_column = "seurat_annotations")
    
    session$setInputs(show_cat_table = T)
    
    session$setInputs(subsample_sce = T,
                      panel_size = 32,
                      display_options = "Marker-marker correlation",
                      select_aa = c("sELISA"),
                      heatmap_expression_norm = "Expression",
                      marker_strategy = "fm")
    
    session$setInputs(start_analysis = T)
    
    expect_false(is.null(cell_types_missing_markers()))
    
  })
})

context("test that setting the user time zone works as intended")

test_that("The user is able to change the time zone properly", {
  testServer(cytosel::cytosel(), expr = {
    
    # expect that the session info is rendered without reactive trigger
    expect_false(is.null(output$current_session_info))
    expect_true(!isTruthy(time_zone_set()))
    
    startup_zone <- Sys.timezone()
    
    # this assumes that the default timezone in the server will never be Kwajalein
    session$setInputs(time_zone_options = "Kwajalein",
                      time_zone = T, pick_time_zone = T)
    
    expect_equal(time_zone_set(), "Kwajalein")
    expect_false(is.null(output$current_time))
    
  })
})

context("test that detection of lowercase genes names registers mouse")

test_that("cytosel is able to find lowercase genes as non-human", {
  testServer(cytosel::cytosel(), expr = {
    
    expect_true(proper_organism())
    
    # read in all lowercase genes for inferorg (should detect as mouse)
    session$setInputs(input_scrnaseq = list(datapath =
                            test_path("pbmc_lowercase.rds")))
    
    expect_false(proper_organism())
    
  })
})

context("test that uploading an RDS in the wrong format generates an error")

test_that("cytosel is able to identify an RDS that is not of the proper SCE format", {
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("fake_rds.rds")))
    expect_false(proceed_with_analysis())
    
  })
})

gc()

# context("test that cytosel can identify multimarkers")
# 
# test_that("cytosel is able to identify multimarkers in a lung dataset 
#           (many cell types for the panel size)", {
#             
#             # skip_on_ci()
#             
#             testServer(cytosel::cytosel(), expr = {
#               
#               session$setInputs(input_scrnaseq = list(datapath =
#                                                         test_path("pbmc_small.rds")),
#                                 coldata_column = "fake_col",
#                                 pick_curated = T,
#                                 min_category_count = 2,
#                                 subset_number = 250)
#               
#               session$setInputs(subsample_sce = F,
#                                 marker_strategy = "fm",
#                                 display_options = "Marker-marker correlation",
#                                 heatmap_expression_norm = "Expression")
#               
#               session$setInputs(tabs = NULL,
#                                 metrics_toggle = NULL,
#                                 select_aa = NULL,
#                                 panel_sorter = "Group by cell type")
#               
#               expect_null(multimarkers())
#               
#               session$setInputs(panel_size = 12,
#                                 start_analysis = T)
#               
#               expect_false(is.null(multimarkers()))
#               
#             })
#           })


test_that("Second Re-upload works on server", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_false(reupload_analysis())
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")))
    
    session$setInputs(panel_size = 24, coldata_column = "seurat_annotations",
                      subsample_sce = F,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression",
                      marker_strategy = "fm")
    
    session$setInputs(read_back_analysis = list(datapath =
                                                  test_path("test_config_2.yml")))
    
    session$setInputs(precomputed_dim = T, select_precomputed_umap = "UMAP",
                      possible_precomputed_dims = reducedDimNames(sce()))
    
    expect_equal(length(specific_cell_types_selected()), 
                 length(unique(sce()[["seurat_annotations"]])))
    
  })
})

context("Check that uploading a minimal yml with just markers works as expected")

#### Re-upload analysis ######
test_that("Re-upload works on server", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_false(reupload_analysis())
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      read_back_analysis = list(datapath =
                                                  test_path("test_config.yml")))
    # verify re-upload analysis reactive worked
    expect_true(reupload_analysis())
    
    # verify that the reactive values were populated from the yml
    expect_equal(length(specific_cell_types_selected()), 5)
    expect_equal(cell_min_threshold(), 10)
    expect_equal(length(markers_reupload()$top_markers), 18)
    expect_equal(length(markers_reupload()$scratch_markers), 6)
    expect_equal(length(specific_cell_types_selected()), 5)
    
    session$setInputs(panel_size = 200, coldata_column = "seurat_annotations",
                      subsample_sce = F,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression",
                      marker_strategy = "fm")
    
    session$setInputs(tabs = NULL,
                      metrics_toggle = NULL,
                      select_aa = NULL,
                      panel_sorter = "Group by cell type")
    
    session$setInputs(start_analysis = T)
    
    expect_equal(length(current_markers()$top_markers),
                 length(markers_reupload()$top_markers))
    
    expect_equal(length(current_markers()$scratch_markers),
                 length(markers_reupload()$scratch_markers))
    
  })
})

gc()

context("Test that resetting the cytosel environment works as expected")


#### Reset analysis ######
test_that("Minimal re-uploads recognizxe markers", {
  
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      coldata_column = "seurat_annotations",
                      read_back_analysis = list(datapath =
                                                  test_path("test_config_3.yml")))
    
    # verify that the reactive values were populated from the yml
    
    expect_equal(length(specific_cell_types_selected()), 
                 length(unique(sce()[["seurat_annotations"]])))
    expect_equal(length(markers_reupload()$top_markers), 18)
    expect_equal(length(markers_reupload()$scratch_markers), 6)
    
  })
})

context("Test that a reupload doesn't occur before the dataset is uploaded")

test_that("Re-uploading won't occur before the dataset is processed", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_false(reupload_analysis())
    
    session$setInputs(read_back_analysis = list(datapath =
                                                  test_path("test_own_panel.txt")))
    
    expect_false(reupload_analysis())
    
  })
})

context("reading a simple list of genes will populate the marker list and search for alias")

test_that("Users may upload a simple list of genes (one per line)", {
  
  markers_attempt <- readLines(test_path("test_own_panel.txt"))
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_false(reupload_analysis())
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      read_back_analysis = list(datapath =
                                                  test_path("test_own_panel.txt")))
    
    expect_true(reupload_analysis())
    expect_true(length(markers_reupload()$top_markers) > 10)
    expect_false(all(markers_attempt %in% markers_reupload()$top_markers))
    expect_true(all(markers_reupload()$top_markers %in% allowed_genes()))
    
    expect_false(all(markers_reupload()$top_markers %in% markers_attempt))
    
  })
})


context("Check that adding or replacing markers with aliases works as intended")

test_that("adding and replacing markers can identify gene aliases", {
  
  marks_to_add <- readLines(test_path("test_own_panel.txt"))
  
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      assay = "logcounts", coldata_column = "seurat_annotations")
    
    session$setInputs(subsample_sce = F,
                      panel_size = 5,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression",
                      marker_strategy = "fm",
                      tabs = NULL,
                      metrics_toggle = NULL,
                      select_aa = NULL,
                      panel_sorter = "Group by cell type",
                      start_analysis = T)
    
    expect_equal(length(current_markers()$top_markers), 5)
    
    session$setInputs(uploadMarkers = list(datapath =
                                             test_path("test_own_panel.txt")),
                      add_to_selected = T)
    
    expect_true(length(current_markers()$top_markers) > 5)
    
    second_marks <- current_markers()$top_markers
    
    expect_false("hCD40L" %in% current_markers()$top_markers)
    expect_true("CD40LG" %in% current_markers()$top_markers)
    
    
    session$setInputs(uploadMarkers = list(datapath =
                                             test_path("test_own_panel.txt")),
                      replace_selected = T)
    
    expect_true(length(second_marks) > length(current_markers()$top_markers))
    
    expect_false("hCD40L" %in% current_markers()$top_markers)
    expect_true("CD40LG" %in% current_markers()$top_markers)
    
    session$setInputs(gene_alias_table_viewer = T)
    
    expect_false(is.null(aliases_table()))
    expect_true(nrow(aliases_table()) > 0)
    
    current_table <- aliases_table()

    
    session$setInputs(uploadMarkers = list(datapath =
                                             test_path("cds_1_to_200.txt")),
                      add_to_selected = T)    
    
    expect_true("CD40LG" %in% aliases_table()$Alias)
    expect_true("hCD40L" %in% aliases_table()$`Original Name`)
    expect_true(nrow(aliases_table()) > nrow(current_table))
    
    current_table_2 <- aliases_table()
    
    session$setInputs(uploadMarkers = list(datapath =
                                             test_path("cds_1_to_200.txt")),
                      replace_selected = T) 
    
    expect_false("hCD40L" %in% aliases_table()$`Original Name`)
    expect_true(nrow(aliases_table()) < nrow(current_table_2))
    
  })
  
})

context("Test basic gene aliases table works")

test_that("A gene alias table is populated when uploading markers", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_true(is.null(aliases_table()))
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")))
    
    session$setInputs(read_back_analysis = list(datapath =
                                                  test_path("cds_1_to_200.txt")))
    
    expect_false(is.null(aliases_table()))
    first_alias <- aliases_table()
    expect_true(is.null(aliases_table_subset()))
    
    session$setInputs(select_aa = c("Protein Array"))
    expect_false(is.null(aliases_table_subset()))
    
    expect_true(nrow(first_alias) > nrow(aliases_table_subset()))
    
    
  })
})


context("Test that starting the help guide does not modify the reactive values")

test_that("Basic help guide functionality", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_null(sce())
    expect_equal(pref_assay(), "logcounts")
    
    session$setInputs(tabs = "inputs")
    session$setInputs(help_guide = T)
    
    expect_null(sce())
    expect_equal(pref_assay(), "logcounts")
    
  })
  
})



