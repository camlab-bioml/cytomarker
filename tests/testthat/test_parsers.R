context("Test that the correct file type is identified")

test_that("passing a path that is not a SCE or Seurat object raises an error", {
  fake_sce <- c(NULL)
  expect_error(read_input_scrnaseq(fake_sce))
})

test_that("cytosel can read in an H%ad file", {
  h5ad_file <- test_path("test_sce.h5ad")
  sce <- read_input_scrnaseq(h5ad_file)
  expect_is(sce, 'SingleCellExperiment')
  expect_equivalent(dim(sce), c(100, 500))
})

context("Testing sce object with Ensembl gene rownames and no rowData")

test_that("Testing sce object with Ensembl gene rownames and no rowData", {
  
  wtc_ensgene_rownames <- readRDS(test_path("wtc-ensgene-rownames_sub.rds"))
  ensembl_hugo_check <- check_rowData_for_hugo(wtc_ensgene_rownames, annotables::grch38)
  expect_is(ensembl_hugo_check, 'list')
  expect_equal(ensembl_hugo_check$status, "no_symbols_found_in_rowdata")
  expect_null(ensembl_hugo_check$genes)
  
  ensembl_row_check <- check_rownames_for_ensembl(wtc_ensgene_rownames,
                                                  annotables::grch38)
  expect_equal(ensembl_row_check$status, "use_ensembl_names")
  expect_false(is.null(ensembl_row_check$sce))
  
  # Do not expect Ensembl in the rownames after converted using annotables
  expect_equal(0, length(rownames(ensembl_row_check)$sce[grepl("ENSG",
                                                               rownames(ensembl_row_check)$sce)]))
  
  # Do not expect the original rownames to be found in annotables
  proportions <- calculate_proportion_in_annotables(
    rownames(wtc_ensgene_rownames), annotables::grch38)
  
  expect_equal(0, proportions$proportion)
  expect_equal(0, proportions$gene_num)
  
  gene_parser <- parse_gene_names(wtc_ensgene_rownames, annotables::grch38)
  expect_is(gene_parser, 'SingleCellExperiment')
  expect_equal(dim(rowData(gene_parser)), c(20935, 0))
  
})

context("Testing sce object with null rownames but gene symbols in rowData")

test_that("Testing sce object with null rownames but gene symbols in rowData", {
  
  wtc_hugo_rowData_genes <- readRDS(test_path("wtc_hugo_rowData_genes_sub.rds"))
  null_row_check <- check_rowData_for_hugo(wtc_hugo_rowData_genes, annotables::grch38)
  expect_is(null_row_check, 'list')
  expect_equal(check_rownames_for_hugo(wtc_hugo_rowData_genes), "rownames_are_null")
  expect_equal(null_row_check$status, "use_rowData_names")
  expect_true("TSPAN6" %in% null_row_check$genes)
  expect_equal(33658, length(null_row_check$genes))
  
  ensembl_row_check <- check_rownames_for_ensembl(wtc_hugo_rowData_genes,
                                                  annotables::grch38)
  expect_is(ensembl_row_check, 'list')
  expect_equal(ensembl_row_check$status, "did_not_find_ensembl")
  expect_null(ensembl_row_check$sce)
  
  proportions <- calculate_proportion_in_annotables(null_row_check$genes,
                                                    annotables::grch38)
  
  expect_true(proportions$proportion > 0.6 & proportions$proportion < 0.7)
  expect_true(proportions$gene_num > 21000 & proportions$gene_num < 21500)
  
  gene_parser <- parse_gene_names(wtc_hugo_rowData_genes, annotables::grch38)
  expect_is(gene_parser, 'SingleCellExperiment')
  expect_equal(dim(rowData(gene_parser)), c(33423, 1))
  
})

context("Testing sce object with null rownames and neither Ensembl or gene symbols in rowData")


test_that("Testing sce object with null rownames and neither Ensembl or gene symbols in rowData", {
  
  wtc_null_rownames <- readRDS(test_path("wtc_null_rownames_sub.rds"))
  null_row_check <- check_rowData_for_hugo(wtc_null_rownames, annotables::grch38)
  expect_is(null_row_check, 'list')
  expect_null(null_row_check$genes)
  expect_equal(check_rownames_for_hugo(wtc_null_rownames), "rownames_are_null")
  
  ensembl_row_check <- check_rownames_for_ensembl(wtc_null_rownames,
                                                  annotables::grch38)
  expect_is(ensembl_row_check, 'list')
  expect_equal(ensembl_row_check$status, "did_not_find_ensembl")
  expect_null(ensembl_row_check$sce)
  
  expect_error(parse_gene_names(wtc_null_rownames, annotables::grch38))
  
})

context("Testing sce object gene symbols in rownames")


test_that("Testing sce object gene symbols in rownames", {
  
  wtc_normal_rownames_factor_celltype <- readRDS(test_path("wtc_normal_rownames_factor_celltype_sub.rds"))
  hugo_row_check <- check_rowData_for_hugo(wtc_normal_rownames_factor_celltype, annotables::grch38)
  expect_is(hugo_row_check, 'list')
  # DO not expect any symbols in the rowData
  expect_null(hugo_row_check$genes)
  expect_equal(hugo_row_check$status, "no_symbols_found_in_rowdata")
  # expect gene symbols in the rownames
  expect_equal(check_rownames_for_hugo(wtc_normal_rownames_factor_celltype,
                                       annotables::grch38),
               "use_rownames")
  
  
  ensembl_row_check <- check_rownames_for_ensembl(wtc_normal_rownames_factor_celltype,
                                                  annotables::grch38)
  # Do not expect to find any human ensembl Ids in the rownames
  expect_is(ensembl_row_check, 'list')
  expect_equal(ensembl_row_check$status, "no_human_ensID")
  expect_null(ensembl_row_check$sce)
  
  gene_parser <- parse_gene_names(wtc_normal_rownames_factor_celltype, annotables::grch38)
  expect_is(gene_parser, 'SingleCellExperiment')
  
  # do not expect any Ensembl Ids in the final parsed sce
  expect_equal(0, length(rownames(gene_parser)[grepl("ENSG",
                                                     rownames(gene_parser))]))
  
})


context("Testing sce object with gene symbol rownames but with a low proportion")

test_that("Testing sce object gene symbols in rownames", {
  
  hugo_rowData_low_proportion <- readRDS(test_path(
    "hugo_rowData_low_proportion.rds"))
  
  rowdata_check <- check_rowData_for_hugo(hugo_rowData_low_proportion,
                                          annotables::grch38)
  
  expect_equal(rowdata_check$status, "low_number_of_gene_matches")
  expect_equal(length(rowdata_check$genes), 33658)
  proper_genes <- rowdata_check$genes[!grepl("Gene", rowdata_check$genes)]
  expect_equal(length(proper_genes), 150)
  
  expect_warning(parse_gene_names(hugo_rowData_low_proportion, annotables::grch38))
  
  rownames(hugo_rowData_low_proportion) <- rowData(hugo_rowData_low_proportion)[,1]
  expect_equal(check_rownames_for_hugo(hugo_rowData_low_proportion, annotables::grch38),
               "low_number_of_gene_matches")
  
  expect_warning(parse_gene_names(hugo_rowData_low_proportion, annotables::grch38))
  expect_equal(dim(parse_gene_names(hugo_rowData_low_proportion, annotables::grch38)),
               c(33657, 200))
  
  rownames(hugo_rowData_low_proportion) <- NULL
  expect_equal(check_rownames_for_hugo(hugo_rowData_low_proportion, annotables::grch38),
               "rownames_are_null")
  expect_warning(parse_gene_names(hugo_rowData_low_proportion, annotables::grch38))
})


context("Testing sce object with Ensembl rownames but with a low proportion")

test_that("Testing sce object Ensembl in rownames", {
  
  ensembl_rownames_low_proportion <- readRDS(test_path(
    "ensembl_rownames_low_proportion.rds"))
  
  expect_equal(0, calculate_proportion_in_annotables(rownames(ensembl_rownames_low_proportion),
                                                     annotables::grch38)$proportion)
  
  rownames_ensembl_check <- check_rownames_for_ensembl(ensembl_rownames_low_proportion,
                                                       annotables::grch38)
  
  expect_equal("low_number_of_gene_matches", rownames_ensembl_check$status)
  expect_equal(dim(rownames_ensembl_check$sce),
               c(101, 200))
  
  gene_parsing <- parse_gene_names(ensembl_rownames_low_proportion, annotables::grch38)
  expect_warning(parse_gene_names(ensembl_rownames_low_proportion, annotables::grch38))
  expect_equal(dim(gene_parsing),
               c(100, 200))
  
  expect_equal(dim(parse_gene_names(ensembl_rownames_low_proportion,
                                    annotables::grch38, remove_confounding_genes = F)),
               c(101, 200))
  
})


context("Testing sce object with partial genes overlaps")

test_that("Partial gene overlaps raise an error", {
  partial_genes <- readRDS(test_path(
    "partial_genes.rds"))
  
  partial_hugo <- check_rownames_for_hugo(partial_genes, annotables::grch38)
  expect_equal(partial_hugo, "did_not_find_hugo_in_rownames")
  expect_error(parse_gene_names(partial_genes, annotables::grch38))
  
})

context("Testing sce object from a different organism")

test_that("Genes from a non-human or mouse organism throws error", {
  diff_organism <- readRDS(test_path(
    "non_human_genes.rds"))
  
  non_ensembl <- check_rownames_for_ensembl(diff_organism, annotables::grch38)
  expect_equal(non_ensembl$status, "did_not_find_ensembl")
  expect_error(parse_gene_names(diff_organism, annotables::grch38))
  
})

context("Test the identification of non-gene elements")

test_that("Ids that have ENS but not ENSG are treated as non-genes", {
  diff_organism <- readRDS(test_path(
    "genes_from_rnor6.rds"))
  
  expect_equal(check_rowData_for_hugo(diff_organism)$status, "no_symbols_found_in_rowdata")
  
  non_ensembl <- check_rownames_for_ensembl(diff_organism, annotables::grch38)
  expect_equal(non_ensembl$status, "no_human_ensID")
  expect_null(non_ensembl$sce)
  expect_error(parse_gene_names(diff_organism, annotables::grch38))
  
})


  
  