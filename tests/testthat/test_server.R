
context("Test basic Shiny app server functionality")

test_that("Server has basic functionality", {
  
  testServer(cytosel::cytosel(), expr = {
    
    # Verify the input sce and default + selected assays
    expect_equal(pref_assay(), "logcounts")
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")))
    
    expect_is(sce(), 'SingleCellExperiment')
    expect_equivalent(dim(sce()), c(13714, 100))
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
                      panel_size = 200, min_category_count = 0,
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
                      add_selected_to_analysis = T,
                      min_category_count = 2,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression",
                      umap_options = "Cell Type",
                      umap_panel_options = NULL,
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
    expect_equivalent(dim(sce()), c(13714, 100))
    
    # do not keep Platelet even though it was not selected for removal
    expect_equivalent(subset(colData(sce()), 
                             seurat_annotations == "Platelet")$keep_for_analysis,
                      "No")
    
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
    
    # expect of 200 markers to get somewhere between 150 to 200 for redundancy
    expect_gt(length(as.character(heatmap()[["x"]][["attrs"]][[1]][["x"]])), 150)
    expect_lte(length(as.character(heatmap()[["x"]][["attrs"]][[1]][["x"]])),
                 200)
    
    # expect that the marker columns are set properly and the lengths are verified
    
    # expect that the marker columns are set properly and the lengths are verified
    session$setInputs(bl_scratch = current_markers()$top_markers[1:10],
                      bl_top = current_markers()$top_markers[11:24])
    
    expect_equal(num_markers_in_selected(), 14)
    expect_equal(num_markers_in_scratch(), 10)
    
    expect_equal(output$scratch_marker_counts, "<B> Scratch Markers: 10 </B>")
    expect_equal(output$selected_marker_counts, "<B> Selected Markers: 14 </B>")
    
    # suggest markers to remove
    
    session$setInputs(suggest_gene_removal = T,
                      n_genes = 3)
    
    expect_equal(length(suggestions()), 3)
    
    # move 5 of the top markers into scratch and assert new numbers
    session$setInputs(markers_to_remove = suggestions(),
                      remove_suggested = T)
    
    # a gene that doesn't exist will not produce suggested replacements
    expect_equal(length(current_markers()$top_markers), 11)
    
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
    
    expect_equal(0, nrow(replacements()[grepl("^RP[L|S]|^MT-|^MALAT",
                                              replacements()$Gene),]))
    
    session$setInputs(alternative_markers_rows_selected = c(1, 3, 4, 7, 8, 10),
                      send_markers = T,
                      send = T)
    
    expect_equal(length(current_markers()$top_markers), 17)
    
    # if you try to add a fake marker, the length will stay the same
    session$setInputs(markers_change_modal = T,
                      add_markers = "FAKE_GENE",
                      enter_marker = T)
    
    expect_equal(num_markers_in_selected(), 17)
    
    
    # if you try to add a marker that's already there, the length will stay the same
    session$setInputs(markers_change_modal = T,
                      add_markers = current_markers()$top_markers[1],
                      enter_marker = T)
    
    expect_equal(num_markers_in_selected(), 17)
    
    # add a gene that isn't in either selected or scratch to increase the count by 1
    different <- allowed_genes()[!allowed_genes() %in% current_markers()$top_markers]
    different <- different[!different %in% current_markers()$scratch_markers]
    
    session$setInputs(markers_change_modal = T,
                      add_markers = different[1],
                      enter_marker = T)
    
    expect_equal(num_markers_in_selected(), 18)
    
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
    
    # after adding 5 markers to 18, assert how many top markers there should be
    expect_equal(length(current_markers()$top_markers), 23)
    
    # if you try to upload a gene that doesn't exist, the length will remain the same
    session$setInputs(uploadMarkers = list(datapath =
                                             test_path("fake_upload.txt")),
                      add_to_selected = T)
    
    # after adding 5 markers to 18, assert how many top markers there should be
    expect_equal(length(current_markers()$top_markers), 23)
    
    
    # Look for the uploaded markers in the top markers
    expect_true("EEF2" %in% current_markers()$top_markers)
    expect_true("MARCKS" %in% current_markers()$top_markers)
    
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

context("Test re-upload and reset Shiny app server functionality")

#### Re-upload analysis ######
test_that("Re-upload and reset works on server", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_false(reupload_analysis())
  
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      read_back_analysis = list(datapath =
                                                  test_path("test_config.yml")))
    # verify reupload analysis reactive worked
    expect_true(reupload_analysis())
    
    # verify that the reactive values were populated from the yml
    expect_equal(length(specific_cell_types_selected()), 6)
    expect_equal(cell_min_threshold(), 10)
    expect_equal(length(markers_reupload()$top_markers), 18)
    expect_equal(length(markers_reupload()$scratch_markers), 6)
    
    session$setInputs(bl_scratch = markers_reupload()$scratch_markers,
                      bl_top = markers_reupload()$top_markers)
    
    session$setInputs(create_reset = T, reset_marker_panel = T)
    
    # resetting the panel sets current markers and input rank lists to 0
    expect_true(reset_panel())
    expect_null(current_markers()$top_markers)
    expect_equal(num_markers_in_selected(), 0)
    expect_equal(num_markers_in_scratch(), 0)
    
  })
})

context("Downloading through the server works as intended")

#### download analysis ######
test_that("Download works on server", {
  
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      read_back_analysis = list(datapath =
                                                  test_path("test_config.yml")))
    
    session$setInputs(panel_size = 24, coldata_column = "seurat_annotations",
                      start_analysis = T)
    
    expect_false("NK" %in% cell_types_to_keep())
    expect_equal(length(cell_types_to_keep()), 5)
    
    expect_false(is.null(output$downloadData))
    
    # test that download feature works with temp dir
    withr::with_tempdir({
      session$setInputs(downloadData = T)
      expect_true(file.exists(output$downloadData))
    })
    
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

context("Test UMAP, violin, and heatmap colouring changes")

test_that("Changing the UMAP, violin, and heatmap colourings work", {
  testServer(cytosel::cytosel(), expr = {
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")),
                      user_selected_cells = c("CD8 T", "Memory CD4 T",
                                              "Naive CD4 T",
                                              "Platelet"),
                      add_selected_to_analysis = T,
                      assay_select = "counts", assay = "counts",
                      coldata_column = "seurat_annotations",
                      min_category_count = 2,
                      subsample_sce = T,
                      panel_size = 20,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression",
                      marker_strategy = "fm")

    session$setInputs(umap_options = "Cell Type", umap_panel_options = "S100A9",
                      start_analysis = T)
    
    heatmap_1 <- heatmap()
    
    # check defaults for UMAP plots
    expect_false(umap_top_gene())
    expect_false(umap_all_gene())
    expect_equal(umap_colouring(), "Cell Type")
    
    session$setInputs(umap_options = "Panel Marker",
    umap_panel_options = "S100A9", umap_panel_cols = T)
    
    # switching the UMAP to gene works
    expect_false(isFALSE(umap_top_gene()))
    expect_false(isFALSE(umap_all_gene()))
    expect_equal(umap_colouring(), "Panel Marker")
    
    viol_markers <- c("EEF2", "RBM3", "MARCKS", "MSN", "JUNB")
    
    # setting the violin plots with genes works
    session$setInputs(genes_for_violin = viol_markers,
                      add_violin_genes = T, viol_viewer = "By Marker")
    expect_false(is.null(output$expression_violin))
    expect_true(all(viol_markers %in% current_markers()$top_markers))
    
    viol_1 <- output$expression_violin
    
    session$setInputs(genes_for_violin = viol_markers,
                      add_violin_genes = T, viol_viewer = "By Cell Type")
    
    expect_false(identical(viol_1, output$expression_violin))
    
    # if set to NULL, still renders an empty violin plot
    session$setInputs(genes_for_violin = NULL)
    expect_false(is.null(output$expression_violin))
    
    # switching the heatmap category makes a different heatmap
    session$setInputs(display_options = "seurat_annotations",
    heatmap_expression_norm = "Expression", start_analysis = T)
    
    expect_false(identical(heatmap_1, heatmap()))
    
  })
})

context("Test the loading of the curated datasets from dropbox")

test_that("Picking the curated dataset works as intended", {
  testServer(cytosel::cytosel(), expr = {
    
    session$setInputs(curated_dataset = T, curated_options = "PBMC small",
                      pick_curated = T)
    expect_true(file.exists("curated/pbmc_small.rds"))
    expect_equivalent(dim(sce()), c(13714, 100))
    dim_1 <- dim(sce())
    
    session$setInputs(curated_dataset = T, curated_options = "PBMC large",
                      pick_curated = T)
    expect_true(file.exists("curated/scRNASeq-test.rds"))
    expect_equivalent(dim(sce()), c(33658, 3220))
    
    expect_false(identical(dim(sce()), dim_1))
    
  })
  
})



