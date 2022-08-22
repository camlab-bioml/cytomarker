
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
    
    
    session$setInputs(user_selected_cells = c("CD8 T", "Memory CD4 T",
                                              "Naive CD4 T",
                                              "Platelet"),
                      panel_size = 24, min_category_count = 0,
                      coldata_column = "seurat_annotations",
                      column = "seurat_annotations",
                      subsample_sce = T,
                      start_analysis = T,
                      add_selected_to_analysis = T,
                      marker_strategy = "fm")
    
    # if the min count is set below 2, do not proceed with analysis
    expect_false(proceed_with_analysis())
    
    expect_null(plots$top_plot)
    
    session$setInputs(min_category_count = 2,
                      start_analysis = T)
    
    # resetting the min count to 2 allows to proceed with analysis
    expect_true(proceed_with_analysis())
    
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
    
    expect_equal(names(cytosel_palette()), c("CD8 T", "Memory CD4 T", "Naive CD4 T"))
    
  })
  
  
})