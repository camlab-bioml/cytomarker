
context("Data reading")

test_that("SingleCellExperiment object reading", {
  obj <- file.path("~/Github/cytosel/tests/testthat/test_sce.rds")
  sce <- read_input_scrnaseq(obj)

  expect_is(sce, 'SingleCellExperiment')
  expect_equivalent(dim(sce), c(100, 500))
  expect_equivalent(rownames(sce), paste("feature-", seq(1, 100, 1), sep=""))
})

test_that("Seurat object reading", {
  obj <- file.path("~/Github/cytosel/tests/testthat/test_seu.rds")
  seu <- read_input_scrnaseq(obj)

  ## Seurat objects are converted to SingleCellExperiments by read_input_scrnaseq
  expect_is(seu, 'SingleCellExperiment')
  expect_equivalent(dim(seu), c(100, 500))
  expect_equivalent(rownames(seu), paste("feature-", seq(1, 100, 1), sep=""))
})

context("Marker finding")

test_that("good_col returns valid columns", {
  sce_path <- system.file("tests/testthat/test_sce.rds", package="cytosel")
  sce <- read_input_scrnaseq(sce_path)
  columns <- good_col(sce, colnames(colData(sce)))

  expect_is(columns, 'list')
  expect_equal(names(columns), c("good", "bad"))
  expect_equal(names(columns$bad), c("colname", "n"))
  expect_equal(length(columns$good), 1)
  expect_equal(columns$good, "col1")
  expect_equal(length(columns$bad$colname), length(columns$bad$n))
  expect_equal(columns$bad$colname, c("orig.ident", "nCount_RNA", "nFeature_RNA", "col2", "ident"))
})

test_that("get_markers and compute_fm returns valid output", {
  sce_path <- system.file("tests/testthat/test_sce.rds", package="cytosel")
  sce <- read_input_scrnaseq(sce_path)
  
  sce$cell_type <- sample(c("A", "B"), ncol(sce), replace=TRUE)
  
  fms <- compute_fm(sce, 
             "cell_type", 
              "logcounts",
            rownames(sce)
  )
  
  expect_is(fms, 'list')
  expect_equal(length(fms), 1)
  expect_equal(nrow(fms[[1]][[1]]), nrow(sce))
  
  markers <- get_markers(fms, panel_size = 32, marker_strategy = 'standard', 
                         sce = sce, 
                         allowed_genes = rownames(sce))


  expect_is(markers, 'list')
  expect_equal(names(markers), c("recommended_markers", "scratch_markers", "top_markers"))
  expect_gt(length(markers$recommended_markers), 0)
  expect_gt(length(markers$top_markers), 0)
})


test_that("get_markers errors for non unique values", {
  obj <- file.path("~/Github/cytosel/tests/testthat/test_sce.rds")
  sce <- read_input_scrnaseq(obj)

  expect_error(get_markers(sce, 'col2', 10, 'logcounts'))
})


context("Plotting")

test_that("get_umap returns valid dataframe", {
  obj <- file.path("~/Github/cytosel/tests/testthat/test_sce.rds")
  sce <- read_input_scrnaseq(obj)
  
  umap_frame_log <- get_umap(sce, "cell_type", "logcounts")
  umap_frame_norm <- get_umap(sce, "cell_type", "counts")
  expect_is(umap_frame_log, 'data.frame')
  expect_equal(ncol(umap_frame_log), 2)
  expect_equal(nrow(umap_frame_log), ncol(sce))
  expect_is(umap_frame_norm, 'data.frame')
  expect_equal(ncol(umap_frame_norm), 2)
  expect_equal(nrow(umap_frame_log), ncol(sce))
  expect_false(unique(umap_frame_norm[,1] == umap_frame_log[,1]))
  expect_false(unique(umap_frame_norm[,2] == umap_frame_log[,2]))

})
# 
# context("Scoring")
# 
# test_that("get_scores returns valid scores", {
# 
# })