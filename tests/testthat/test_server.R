
if (!interactive()) {
  skip("Do not run server side tests if non-interactive")
  
}

skip_on_cran()

context("Test Shiny app server functionality")

test_that("Server has functionality", {
  
  # appfile <- system.file("app.R", package = "cytosel")
  
  testServer(cytosel::cytosel(), expr = {
    
    expect_equal(pref_assay(), "logcounts")
    session$setInputs(input_scrnaseq = list(datapath =
                                              system.file("/tests/testthat/test_sce_with_real_fake_gene_symbol.rds",
                                                          package="cytosel",
                                                          lib.loc = "cytosel",
                                                          mustWork = T)))
    expect_is(sce(), 'SingleCellExperiment')
    expect_equivalent(dim(sce()), c(100, 500))
    session$setInputs(assay_select = "counts", assay = "counts")
    expect_equal(pref_assay(), "counts")
    expect_equal(output$selected_assay, paste("Selected assay: ", "counts"))
  })
  
  
})