
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
    
    session$setInputs(user_selected_cells = c("CD8 T", "Memory CD4 T",
                                                       "Naive CD4 T",
                                              "Platelet"),
                      panel_size = 24, min_category_count = 2,
                      coldata_column = "seurat_annotations",
                      proceed_with_analysis = T,
                      start_analysis = T)
    
    expect_equivalent(unique(sce()$keep_for_analysis), c("No", "Yes"))
    expect_equivalent(dim(sce()), c(13714, 100))
    
    # do not keep Platelet even though it was selected as there is only one count
    expect_equivalent(subset(colData(sce()), 
                             seurat_annotations == "Platelet")$keep_for_analysis,
                      "No")
    
  })
  
  
})