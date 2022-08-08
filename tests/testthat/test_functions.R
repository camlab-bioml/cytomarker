
context("Data reading")

# test_that("SingleCellExperiment object reading", {
#   obj <- readRDS("~/Github/cytosel/tests/testthat/test_sce.rds")
#   sce <- read_input_scrnaseq(obj)
# 
#   expect_is(sce, 'SingleCellExperiment')
# })
# 
# test_that("Seurat object reading", {
#   obj <- readRDS("~/Github/cytosel/tests/testthat/test_seu.rds")
#   seu <- read_input_scrnaseq(obj)
# 
#   ## Seurat objects are converted to SingleCellExperiments by read_input_scrnaseq
#   expect_is(seu, 'SingleCellExperiment')
# })

context("Marker finding")

test_that("good_col returns valid columns", {
  sce_path <- system.file("tests/testthat/test_sce.rds", package="cytosel")
  sce_path <- "test_sce.rds"
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
  sce_path <- "test_sce.rds"
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
# 
# 
# test_that("get_markers errors for non unique values", {
#   obj <- readRDS("~/Github/cytosel/tests/testthat/test_sce.rds")
#   sce <- read_input_scrnaseq(obj)
#   
#   expect_error(get_markers(sce, 'col2', 10, 'logcounts'))
# })

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

context("Testing palette")

test_that("palette always returns consistent colors with and without seeding", {
  test_palettes <- list()
  
  pal_1 = c("#7FC97F", "#FFFF99", "#386CB0", "#6A3D9A","#CCEBC5")
  pal_2 <- c("#BEAED4", "#666666", "#A6CEE3", "#CAB2D6", "#B3CDE3",
             "#66C2A5", "#A6D854")
  
  pal_seed_1 <- c("#CCEBC5", "#80B1D3", "#E6F5C9", "#FDCDAC", "#CBD5E8")
  pal_seed_2 <- c("#FFFF33", "#BC80BD", "#6A3D9A", "#FDB462", "#FC8D62",
                  "#E31A1C", "#FFED6F")
  
  i <- 1
  while(i <= 10) {
    new_pal <- create_global_colour_palette()
    new_pal_1 <- as.vector(new_pal[c(1, 4, 5, 26, 73)])
    new_pal_2 <- as.vector(new_pal[c(2, 8, 17, 25, 30, 55, 59)])
    
    
    expect_equal(new_pal_1, pal_1)
    expect_equal(new_pal_2, pal_2)
    
    
    new_pal_seed <- create_global_colour_palette(pal_seed = 123)
    new_pal_1_seed <- as.vector(new_pal_seed[c(1, 4, 5, 26, 73)])
    new_pal_2_seed <- as.vector(new_pal_seed[c(2, 8, 17, 25, 30, 55, 59)])
    
    expect_equal(new_pal_1_seed, pal_seed_1)
    expect_equal(new_pal_2_seed, pal_seed_2)
    
    i <- i + 1
  }
  
})
