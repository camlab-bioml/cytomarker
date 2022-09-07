library(yaml)
context("Test basic functions")

test_that("round3 function is effective", {
  expect_equal(round3(0.3456), "0.300")
  obj <- test_path("pbmc_small.rds")
  sce <- read_input_scrnaseq(obj)
  expect_is(convert_column_to_character_or_factor(sce, "num_col")$num_col,
            "character")
  
})

context("Data reading")

test_that("SingleCellExperiment object reading", {
  obj <- test_path("test_sce.rds")
  sce <- read_input_scrnaseq(obj)

  expect_is(sce, 'SingleCellExperiment')
  expect_equivalent(dim(sce), c(100, 500))
  expect_equivalent(rownames(sce), paste("feature-", seq(1, 100, 1), sep=""))
})

test_that("Seurat object reading", {
  obj <- test_path("test_seu.rds")
  seu <- read_input_scrnaseq(obj)

  ## Seurat objects are converted to SingleCellExperiments by read_input_scrnaseq
  expect_is(seu, 'SingleCellExperiment')
  expect_equivalent(dim(seu), c(100, 500))
  expect_equivalent(rownames(seu), paste("feature-", seq(1, 100, 1), sep=""))
})

context("Marker finding")

test_that("good_col returns valid columns", {
  sce_path <- test_path("test_sce.rds")
  sce <- read_input_scrnaseq(sce_path)
  columns <- good_col(sce, colnames(colData(sce)))

  expect_is(columns, 'list')
  expect_equal(names(columns), c("good", "bad"))
  expect_equal(names(columns$bad), c("colname", "n"))
  expect_equal(length(columns$good), 1)
  expect_equal(columns$good, "col1")
  expect_equal(length(columns$bad$colname), length(columns$bad$n))
  expect_equal(columns$bad$colname, 
               c("orig.ident", "nCount_RNA", "nFeature_RNA", "col2", "ident"))
})

test_that("get_markers and compute_fm returns valid output", {
  sce_path <- test_path("test_sce.rds")
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
  
  
  markers_geneBasis <- get_markers(fms, panel_size = 24, marker_strategy = 'geneBasis',
                                   sce = sce,
                                   allowed_genes = rownames(sce))
  
  expect_is(markers, 'list')
  expect_equal(names(markers_geneBasis), 
               c("recommended_markers", "scratch_markers", "top_markers"))
  expect_null(markers_geneBasis$scratch_markers)
  expect_gt(length(markers_geneBasis$recommended_markers), 0)
  expect_gt(length(markers_geneBasis$top_markers), 0)
  
  # expect the different selections to give different top markers
  expect_false(setequal(sort(markers_geneBasis$topmarkers),
                        sort(markers$top_markers)))
  
})

test_that("get_markers errors for non unique values", {
  obj <- test_path("test_sce.rds")
  sce <- read_input_scrnaseq(obj)

  expect_error(get_markers(sce, 'col2', 10, 'logcounts'))
})


context("Plotting")

test_that("get_umap returns valid dataframe and values with different assays", {
  obj <- test_path("pbmc_small.rds")
  sce <- read_input_scrnaseq(obj)

  umap_frame_log <- get_umap(sce, "seurat_annotations", "logcounts")
  umap_frame_norm <- get_umap(sce, "seurat_annotations", "counts")
  expect_is(umap_frame_log, 'data.frame')
  expect_equal(ncol(umap_frame_log), 3)
  expect_equal(nrow(umap_frame_log), ncol(sce))
  expect_is(umap_frame_norm, 'data.frame')
  expect_equal(ncol(umap_frame_norm), 3)
  expect_equal(nrow(umap_frame_log), ncol(sce))
  expect_false(unique(umap_frame_norm[,1] == umap_frame_log[,1]))
  expect_false(unique(umap_frame_norm[,2] == umap_frame_log[,2]))

})

context("heatmaps")

test_that("create_heatmap works effectively with different normalizations", {
  obj <- test_path("pbmc_small.rds")
  sce <- read_input_scrnaseq(obj)
  
  sce$num_placeholder <- rep(seq(5), 20)

  markers <- list(top_markers = c("EEF2", "RBM3", "MARCKS", "MSN", "JUNB"))
  
  heat_z <- create_heatmap(sce, markers, "num_placeholder", "Expression",
                 "z-score", "logcounts")
  
  expect_is(heat_z, 'plotly')
  expect_equal(as.character(heat_z$x$attrs[[1]]$z)[2], "expression")
  
  heat_cor <- create_heatmap(sce, markers, "seurat_annotations", "Marker-marker correlation",
                           "z-score", "logcounts")
  
  expect_is(heat_cor, 'plotly')
  expect_false(as.character(heat_cor$x$attrs[[1]]$z)[2] == "expression")
  
})

# context("Scoring")
#
# test_that("get_scores returns valid scores", {
#
# })
  

#   
# })

context("Testing palette and colour conversions")

test_that("Palette always returns consistent colors with and without seeding", {
  test_palettes <- list()
  
  pal_1 = c("#88CCEE", "#117733", "#332288", "#E6AB02", "#E5C494")
  pal_2 <- c("#CC6677", "#999933", "#386CB0", "#66A61E",
             "#1F78B4", "#FFF2AE","#377EB8")
  
  pal_seed_1 <- c("#88CCEE", "#117733", "#332288", "#FBB4AE", "#33A02C")
  pal_seed_2 <- c("#CC6677", "#999933", "#E6F5C9", "#1B9E77",
                  "#BF5B17", "#386CB0", "#1F78B4")
  
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

test_that("Conversion of the background colours always yields black or white",
          {
            background_texts <- sapply(create_global_colour_palette(), 
            FUN = function(x) set_text_colour_based_on_background(x))
            
            expect_equal(length(background_texts[!is.null(background_texts)]), 100)
            expect_equivalent(unique(background_texts), c("#000000", "#ffffff"))
            
          })


context("Check filtering steps")

test_that("Filtering sce objects by minimum count passes & retains original size", {
  
  obj <- test_path("pbmc_small.rds")
  sce <- read_input_scrnaseq(obj)
  
  grouped <- create_table_of_hetero_cat(
    sce, "seurat_annotations")
  expect_is(grouped, 'data.frame')
  expect_equal(nrow(grouped), 9)
  
  cells_retained <- remove_cell_types_by_min_counts(
  grouped, sce, "seurat_annotations", 5)
  
  expect_vector(cells_retained)
  
  # expect 3 cell types dropped
  expect_equal(length(cells_retained), 6)
  # expect cells types below 5 are not retained
  expect_false("Platelet" %in% cells_retained)
  expect_false("Undetermined" %in% cells_retained)
  expect_false("FCGR3A+ Mono" %in% cells_retained)
  expect_true("NK" %in% cells_retained)
  
  
  sce_with_retention <- create_sce_column_for_analysis(sce, cells_retained, 
                                                       "seurat_annotations")
  
  expect_is(sce_with_retention, 'SingleCellExperiment')
  # expect that the length of the sce does not change when adding the keep_for_analysis
  # annotations
  expect_equal(dim(sce_with_retention)[2], dim(sce)[2])
  
  # Expect 7 cells not kept for analysis due to filtering
  expect_equivalent(c(7, 93),
                    as.data.frame(table(sce_with_retention$keep_for_analysis))[,2])
  
})


context("Antibody finding in the abcam table")

test_that("Filtering sce objects by minimum count passes & retains original size", {
  antibody_info <- dplyr::rename(cytosel_data$antibody_info, Symbol = `Gene Name (Upper)`)
  antibody_info <- tidyr::drop_na(antibody_info)
  
  
  expect_equal(get_antibody_info("CD45", antibody_info)$status, "NO_GENE_FOUND")
  expect_equal(get_antibody_info("CD74", antibody_info)$status, "ANTIBODY_FOUND")
  
  
})


context("Return types from catch-all error/notification function")

test_that("throw_error_or_warning returns the correct type", {
  antibody_info <- dplyr::rename(cytosel_data$antibody_info, Symbol = `Gene Name (Upper)`)
  antibody_info <- tidyr::drop_na(antibody_info)
  
  expect_error(throw_error_or_warning(type = 'error', message = "Testing error"))
  expect_error(throw_error_or_warning(type = 'notification', 
                                                   message = "Testing notif"))
})

context("test download")

test_that("download works as expected", {
  
  td <- tempdir(check = TRUE)
  obj <- test_path("pbmc_small.rds")
  sce <- read_input_scrnaseq(obj)
  
  heatmap <- create_heatmap(sce, list(top_markers = rownames(sce)[1:100]),
                            "seurat_annotations", "Marker-marker correlation",
                            "z-score", "logcounts")
  
  umap_all <- get_umap(sce, "seurat_annotations", "logcounts")
  umap_top <- get_umap(sce[rownames(sce)[1:100]], "seurat_annotations", "logcounts")
  
  plots <- list()
  
  plots$all_plot <- plot_ly(umap_all, x=~UMAP1, y=~UMAP2,
                            type='scatter', hoverinfo="text") %>% 
    layout(title = "UMAP all genes")
  
  plots$top_plot <- plot_ly(umap_top, x=~UMAP1, y=~UMAP2, 
                            type='scatter', hoverinfo="text") %>% 
    layout(title = "UMAP selected markers")
  
  score_frame <- data.frame(what = c("sam_1", "sam_2"),
                            score = c(0.55, 0.66))
  
  
  plots$metric_plot <- plot_ly(score_frame, x = ~score, y = ~what, type='box', hoverinfo = 'none') %>% 
    layout(xaxis = list(title="Score"),
           yaxis = list(title="Source"))
  
  withr::with_tempdir({
    filepath <- file.path(paste0("Cytosel-Panel-", Sys.Date(), ".zip"))
    
    download_data(filepath,
                  list(top_markers = rownames(sce)[1:100]), 
                  plots, heatmap, "fake_path_to_sce", "logcounts",
                  "seurat_annotations", 24, 2, "no")
    
    # unzip to tempdir and read back
    unzip(filepath, exdir = td)
    
    marks_back <- read.table(file.path(td, paste0("markers-", Sys.Date(), ".txt")))
    expect_equal(marks_back$V1, rownames(sce)[1:100])
    
    yaml_back <- read_yaml(file.path(td, paste0("config-", Sys.Date(), ".yml")))
    expect_equal(yaml_back$`Input file`, "fake_path_to_sce")
    expect_equal(yaml_back$`Heterogeneity source`, "seurat_annotations")
    
  })
  
})

context("test error warnings from modals")

test_that("Error modals throw errors", {
  expect_error(invalid_modal())
  fake_long_cols <- list(colnames = "Too_Many", n = 200)
  expect_error(unique_element_modal(fake_long_cols))
  fake_one_col <- list(colnames = "CAMLAB", n = 1)
  expect_error(unique_element_modal(fake_one_col))
  
  expect_error(warning_modal("CAMLAB"))
  expect_error(dne_modal("CAMLAB"))
  
  expect_error(too_large_to_show_modal("FAKE_COL"))
  
  expect_error(throw_error_or_warning(message = "Error!",
                                        notificationType = 'error'))
  
})










