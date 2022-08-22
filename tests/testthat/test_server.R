# # 
# if (!interactive()) {
#   skip("Do not run server side tests if non-interactive")
# 
# }

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
                      panel_size = 24, min_category_count = 2,
                      coldata_column = "seurat_annotations",
                      column = "seurat_annotations",
                      subsample_sce = T,
                      proceed_with_analysis = T,
                      start_analysis = T)
    
    expect_true("Yes" %in% unique(sce()$keep_for_analysis))
    expect_true("No" %in% unique(sce()$keep_for_analysis))
    expect_equal(length(unique(sce()$keep_for_analysis)), 2)
    expect_equivalent(dim(sce()), c(13714, 100))
    
    # do not keep Platelet even though it was not selected for removal
    expect_equivalent(subset(colData(sce()), 
                             seurat_annotations == "Platelet")$keep_for_analysis,
                      "No")
    
  })
  
  
})