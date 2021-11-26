
context("Data reading")

test_that("SingleCellExperiment object reading", {
  sce <- read_input_scrnaseq("test_sce.rds")
  
  expect_is(sce, 'SingleCellExperiment')
})

test_that("Seurat object reading", {
  seu <- read_input_scrnaseq("test_seu.rds")
  
  ## Seurat objects are converted to SingleCellExperiments by read_input_scrnaseq
  expect_is(seu, 'SingleCellExperiment')
})

context("Marker finding")

# test_that("get_markers returns valid input", {
#   sce <- read_input_scrnaseq("test_sce.rds")
#   markers <- get_markers(sce, 'col1', 10, 'logcounts')
#   
#   expect_is(markers, 'list')
#   expect_equal(names(markers), c("recommended_markers", "scratch_markers", "top_markers"))
#   expect_gt(length(markers$recommended_markers), 0)
#   expect_gt(length(markers$top_markers), 0)
# })


test_that("get_markers errors for non unique values", {
  sce <- read_input_scrnaseq("test_sce.rds")
  
  expect_error(get_markers(sce, 'col2', 10, 'logcounts'))
})

# 
# context("Plotting")
# 
# test_that("get_umap returns valid dataframe", {
#   
# })
# 
# context("Scoring")
# 
# test_that("get_scores returns valid scores", {
#   
# })