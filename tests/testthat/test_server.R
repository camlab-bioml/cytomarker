
if (!interactive()) {
  skip("Do not run server side tests if non-interactive")

}

skip_on_cran()

context("Test Shiny app server functionality")

test_that("Server has functionality", {
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_equal(pref_assay(), "logcounts")
    session$setInputs(input_scrnaseq = list(datapath =
                                              system.file("/tests/testthat/pbmc_small.rds",
                                                          package="cytosel",
                                                          lib.loc = "cytosel",
                                                          mustWork = T)))
    expect_is(sce(), 'SingleCellExperiment')
    expect_equivalent(dim(sce()), c(13714, 100))
    session$setInputs(assay_select = "counts", assay = "counts")
    expect_equal(pref_assay(), "counts")
    expect_equal(output$selected_assay, paste("Selected assay: ", "counts"))
    
    session$setInputs(user_removed_cells = c("Undetermined"),
                      panel_size = 24, min_category_count = 2,
                      coldata_column = "seurat_annotations",
                      column = "seurat_annotations",
                      subsample_sce = T,
                      proceed_with_analysis = T,
                      start_analysis = T)
    
    expect_true("Yes" %in% unique(sce()$keep_for_analysis))
    expect_true("No" %in% unique(sce()$keep_for_analysis))
    expect_equivalent(dim(sce()), c(13714, 100))
    
    # do not keep Platelet even though it was not selected for removal
    expect_equivalent(subset(colData(sce()), 
                             seurat_annotations == "Platelet")$keep_for_analysis,
                      "No")
    
  })
  
  
})