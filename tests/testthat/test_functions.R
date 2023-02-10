library(yaml)
context("Test basic functions")

test_that("conversion functions are effective", {
  expect_equal(round3(0.3456), "0.300")
  obj <- test_path("pbmc_small.rds")
  sce <- read_input_scrnaseq(obj)
  expect_is(convert_column_to_character_or_factor(sce, "num_col")$num_col,
            "character")
  
  expect_true(check_for_human_genes(sce))
  obj <- test_path("pbmc_lowercase.rds")
  sce <- read_input_scrnaseq(obj)
  
  expect_false(check_for_human_genes(sce))
})

context("Pre-processing the antibody applications")

test_that("Processing the antibody applications produces the intended data structures", {
  
  applications_parsed <- get_antibody_applications(cytosel_data$antibody_info, 
                                                   'Symbol', 'Listed Applications')
  
  expect_is(applications_parsed, 'list')
  expect_equal(names(applications_parsed), c("unique_applications", "application_gene_map"))
  expect_equal(length(applications_parsed$unique_applications), 12)
  expect_equal(length(applications_parsed$application_gene_map), 12)
  
})


context("Data reading")

test_that("SingleCellExperiment object reading", {
  obj <- test_path("test_sce.rds")
  sce <- read_input_scrnaseq(obj)
  
  expect_is(sce, 'SingleCellExperiment')
  expect_equivalent(dim(sce), c(100, 500))
  expect_equivalent(rownames(sce), paste("feature-", seq(1, 100, 1), sep=""))
  
  expect_null(read_input_scrnaseq(test_path("fake_rds.rds")))
  
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
  sce_path <- test_path("pbmc_small.rds")
  sce <- read_input_scrnaseq(sce_path)

  # sce$cell_type <- sample(c("A", "B"), ncol(sce), replace=TRUE)

  fms <- compute_fm(sce,
             "seurat_annotations",
              "logcounts",
            rownames(sce)
  )
  
  expect_is(fms, 'list')
  expect_equal(length(fms), 1)
  expect_equal(nrow(fms[[1]][[1]]), nrow(sce))

  markers <- get_markers(fms, panel_size = 100, marker_strategy = 'standard',
                         sce = sce,
                         allowed_genes = rownames(sce))
  
  markers <- markers$marker

  expect_is(markers, 'list')
  expect_equal(names(markers), c("recommended_markers", "scratch_markers", "top_markers"))
  expect_gt(length(markers$recommended_markers), 0)
  expect_gt(length(markers$top_markers), 0)
  
  
  markers_geneBasis <- get_markers(fms, panel_size = 24, marker_strategy = 'geneBasis',
                                   sce = sce,
                                   allowed_genes = rownames(sce))
  
  markers_geneBasis <- markers_geneBasis$marker
  
  expect_is(markers_geneBasis, 'list')
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

  umap_frame_log <- get_umap(sce, "seurat_annotations", "logcounts", 
                             markers_to_use = rownames(sce))
  umap_frame_norm <- get_umap(sce, "seurat_annotations", "counts", 
                              markers_to_use = rownames(sce))
  expect_is(umap_frame_log, 'data.frame')
  expect_equal(ncol(umap_frame_log), 3)
  expect_equal(nrow(umap_frame_log), ncol(sce))
  expect_is(umap_frame_norm, 'data.frame')
  expect_equal(ncol(umap_frame_norm), 3)
  expect_equal(nrow(umap_frame_log), ncol(sce))
  expect_false(unique(umap_frame_norm[,1] == umap_frame_log[,1]))
  expect_false(unique(umap_frame_norm[,2] == umap_frame_log[,2]))
  
  expect_equal(detect_umap_dims_in_sce(sce), "UMAP")

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
  
  heat_z_expression <- heat_z$x$attrs[[3]]$z
  heat_z_expression <- round(heat_z_expression, 3)
  
  expect_equal(heat_z_expression,
               structure(c(-1.049, -1.079, 0.33, 0.839, 0.959, 0.483, -1.289, 
                           0.161, 1.289, -0.645, 0.943, -1.414, -0.471, 0, 0.943, 0.153, 
                           -0.614, -0.614, -0.614, 1.687, 1.356, 0.574, -0.469, -1.252, 
                           -0.209), .Dim = c(5L, 5L), .Dimnames = list(c("3", "1", "5", 
                                                                         "4", "2"), 
                                                                       c("JUNB", "EEF2", 
                                                                         "MSN", "MARCKS", "RBM3"))))
  
  heat_cor <- create_heatmap(sce, markers, "seurat_annotations", "Marker-marker correlation",
                           "z-score", "logcounts")
  
  expect_is(heat_cor, 'plotly')
  expect_false(as.character(heat_cor$x$attrs[[1]]$z)[2] == "expression")
  
})

context("violin plotting") 

test_that("Violin plotting returns the appropriate data structure", {
  obj <- test_path("pbmc_small.rds")
  sce <- read_input_scrnaseq(obj)
  
  markers <- c("EEF2", "RBM3", "MARCKS", "MSN", "JUNB")
  
  viol <- make_violin_plot(sce, markers, "seurat_annotations", "logcounts")
  expect_is(viol, 'data.frame')
  expect_equal(ncol(viol), 4)
  expect_equal(length(unique(viol$Gene)), 5)
  expect_equal(colnames(viol), c("Gene", "Expression",
                                 "Cell Type", "Colour"))
  
})



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
  
  
  # expect that the length of the sce does not change when adding the keep_for_analysis
  # annotations
  sce_no_null <- remove_null_and_va_from_cell_cat(sce, "seurat_annotations")
  
  expect_equal(dim(sce_no_null)[2], dim(sce)[2])
  
  sce_with_retention <- create_sce_column_for_analysis(sce_no_null, cells_retained, 
                                                       "seurat_annotations")
  
  expect_is(sce_with_retention, 'SingleCellExperiment')
  expect_equal(dim(sce_with_retention)[2], dim(sce)[2])
  
  # Expect 7 cells not kept for analysis due to filtering
  expect_equivalent(c(7, 93),
                    as.data.frame(table(sce_with_retention$keep_for_analysis))[,2])
  
})


context("Antibody finding in the abcam table")

test_that("Filtering sce objects by minimum count passes & retains original size", {
  antibody_info <- cytosel_data$antibody_info
  
  
  expect_equal(get_antibody_info("CD45", antibody_info)$status, "NO_GENE_FOUND")
  expect_equal(get_antibody_info("CD74", antibody_info)$status, "ANTIBODY_FOUND")
  
  
})


context("Return types from catch-all error/notification function")

test_that("throw_error_or_warning returns the correct type", {
  antibody_info <- cytosel_data$antibody_info
  
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
  
  umap_all <- get_umap(sce, "seurat_annotations", "logcounts", markers_to_use = rownames(sce))
  umap_top <- get_umap(sce[rownames(sce)[1:100]], "seurat_annotations", "logcounts",
                       markers_to_use = rownames(sce))
  
  plots <- list()
  
  plots$all_plot <- plot_ly(umap_all, x=~UMAP_1, y=~UMAP_2,
                            type='scatter', hoverinfo="text") %>% 
    layout(title = "UMAP all genes")
  
  plots$top_plot <- plot_ly(umap_top, x=~UMAP_1, y=~UMAP_2, 
                            type='scatter', hoverinfo="text") %>% 
    layout(title = "UMAP selected markers")
  
  score_frame <- data.frame(what = c("sam_1", "sam_2"),
                            score = c(0.55, 0.66))
  
  
  plots$metric_plot <- plot_ly(score_frame, x = ~score, y = ~what, type='box', hoverinfo = 'none') %>% 
    layout(xaxis = list(title="Score"),
           yaxis = list(title="Source"))
  
  withr::with_tempdir({
    filepath <- file.path(paste0("Cytosel-Panel-", Sys.Date(), ".zip"))
    
    placeholder_markers <- c("EEF2", "RBM3", "MARCKS", "MSN", "FTL")
    
    fake_table <- cytosel_data$antibody_info |> dplyr::filter(Symbol %in% 
                              placeholder_markers) |>
      mutate(`Host Species` = factor(`Host Species`),
             `Product Category Tier 3` = factor(`Product Category Tier 3`),
             `KO Status` = factor(`KO Status`),
             `Clone Number` = factor(`Clone Number`),
             `Human Protein Atlas` = "fake_link",
             `External Link` = paste0('<a href="',`Datasheet URL`, '"',
                                      ' target="_blank" rel="noopener noreferrer"',
                                      '>',"View in Abcam website",'</a>')) |>
      dplyr::select(-c(`Datasheet URL`))

    markdown_report_path <- system.file(file.path("report", "rmarkdown-report.Rmd"), 
                                        package = "cytosel")
    
    fake_metrics <- data.frame(`Cell Type` = c("Fake_1", "Fake_2", "Fake_3"),
                               Score = c(0.99, 1, 0.8))
    
    base_config <- create_run_param_list(marker_list = list(top_markers = rownames(sce)[1:100]), 
                                         "fake_path_to_sce", "logcounts",
                                         "seurat_annotations", 24, 2, "no",
                                         80,
                                         "fm", NULL, NULL, FALSE, 100, 13714, fake_metrics)
    
    expect_is(base_config, 'list')
    
    fake_genes <- c("Fake_1", "Fake_2", "Fake_3")
    names(fake_genes) <- c("Gene_1", "Gene_2", "gene_3")
    
    fake_overall <- list(score = 0.99, counts = 100)
    
    skip_on_ci()
    
    download_data(filepath, base_config, plots, heatmap, fake_table, markdown_report_path,
                  fake_metrics, fake_overall, fake_genes)
    
    # unzip to tempdir and read back
    unzip(filepath, exdir = td)
    
    yaml_back <- read_yaml(file.path(td, paste0("config-", Sys.Date(), ".yml")))
    expect_equal(yaml_back$`Input file`, "fake_path_to_sce")
    expect_equal(yaml_back$`Heterogeneity source`, "seurat_annotations")
    
  })
  
})

context("test parsing multimarkers")

test_that("function for parsing multimarkers works as intended", {
  
  # skip_on_ci()
  initial <- read.table(test_path("recommendations.tsv"), header = T, sep = "\t")
  recommended <- read.table(test_path("recommended.txt"))$V1
  multimarkers <- create_multimarker_frame(initial, recommended)
  expect_equal(nrow(multimarkers), 1)
  expect_equal(multimarkers$marker, "S100A6")
  
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
  
  expect_null(high_cell_number_warning(99, 100))
  expect_error(high_cell_number_warning(200, 100))
  
  expect_error(reupload_failed_modal())
  expect_error(reupload_before_sce_modal())
  expect_error(reupload_warning_modal("title","body"))
  expect_error(current_pan_not_valid_modal("GENE"))
  
  # expected class from a modal dialog box
  expect_is(reset_analysis_modal(), 'shiny.tag')
  expect_is(suggestion_modal(failed = T, c("Sug_1", "Sug_2"), "Sug_1"), 'shiny.tag')
  expect_is(curated_dataset_modal(c("cur_1", "cur_2"), c("cur_1", "cur_2"), 
                                  c("cur_1", "cur_2"), failed = T), 'shiny.tag')
  expect_error(subsampling_error_modal(c("Type_1", "Type_2")))
  expect_is(time_zone_modal(cytosel_data$time_zones, NULL), 'shiny.tag')
  expect_is(reset_option_on_change_modal("placeholder"), 'shiny.tag')
  expect_error(invalid_modal())
  expect_error(invalid_metadata_modal("fake_column"))
})


