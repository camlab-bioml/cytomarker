
context("Test Shiny app server functionality")

test_that("Server has functionality", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_equal(pref_assay(), "logcounts")
    session$setInputs(input_scrnaseq = list(datapath =
                                              test_path("pbmc_small.rds")))
    expect_is(sce(), 'SingleCellExperiment')
    expect_equivalent(dim(sce()), c(13714, 100))
    session$setInputs(assay_select = "counts", assay = "counts")
    expect_equal(pref_assay(), "counts")
    expect_equal(output$selected_assay, paste("Selected assay: ", "counts"))
    
    session$setInputs(coldata_column = "seurat_annotations",
                      min_category_count = 2,
                      show_cat_table = T)
    
    expect_equal(nrow(summary_cat_tally()), 9)
    expect_equal(length(specific_cell_types_selected()), 9)
    expect_equal(specific_cell_types_selected(), unique(sce()[["seurat_annotations"]]))                
    
    session$setInputs(user_selected_cells = c("CD8 T", "Memory CD4 T",
                                              "Naive CD4 T",
                                              "Platelet"),
                      panel_size = 24, min_category_count = 0,
                      subsample_sce = T,
                      start_analysis = T,
                      add_selected_to_analysis = T,
                      marker_strategy = "fm")
    
    # if the min count is set below 2, do not proceed with analysis
    expect_false(proceed_with_analysis())
    
    expect_null(plots$top_plot)
    
    session$setInputs(min_category_count = 2,
                      display_options = "Marker-marker correlation",
                      heatmap_expression_norm = "Expression",
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
    
    expect_equivalent(dim(sce()), c(13714, 100))
    
    # do not keep Platelet even though it was not selected for removal
    expect_equivalent(subset(colData(sce()), 
                             seurat_annotations == "Platelet")$keep_for_analysis,
                      "No")
    
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
    
    # espect that the marker columns are set properly and the lengths are verified
    session$setInputs(bl_scratch = current_markers()$top_markers[1:10],
                      bl_top = current_markers()$top_markers[11:24],
                      refresh_marker_counts = T)
    
    expect_equal(num_markers_in_selected(), 14)
    expect_equal(num_markers_in_scratch(), 10)
    
    # if you try to add a marker that's already there, the length will stay the same
    session$setInputs(add_markers = current_markers()$top_markers[1],
                      enter_marker = T)
    expect_equal(num_markers_in_selected(), 14)
    
    # add a gene that isn't in either selected or scratch to increase the count by 1
    different <- allowed_genes()[!allowed_genes() %in% current_markers()$top_markers]
    different <- different[!different %in% current_markers()$scratch_markers]
    
    session$setInputs(add_markers = different[1],
                      enter_marker = T)
    expect_equal(num_markers_in_selected(), 15)
    
    # expect null since no markers have been suggested
    expect_true(is.null(selected_cell_type_markers()))
    
    # generate suggested replacement markers
    session$setInputs(input_gene = current_markers()$top_markers[3],
                      number_correlations = 10,
                      enter_gene = T)
    
    # replacements should be of length 10
    expect_false(is.null(replacements()))
    
    expect_equal(nrow(replacements()), 10)
    
    # get 50 suggested markers for the specific cell type
    session$setInputs(cell_type_markers = "CD8 T",
                      add_cell_type_markers = T)
    
    expect_equal(nrow(marker_suggestions()), 50)
    
   expect_false(is.null(output$cell_type_marker_reactable))
   
    
   session$setInputs(add_select_markers = T)
    
   session$setInputs(min_category_count = 50,
                      start_analysis = T)
    
  expect_false(any_cells_present())
    
  })
})