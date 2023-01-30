library(SingleCellExperiment)
library(Matrix)
# library(tidyverse)
library(scater)
library(scran)

set.seed(23523L)

counts <- readMM("data/counts.umi.txt.gz")
cells <- read.table("data/cells.umi.new.txt")
genes <- read.table("data/genes.umi.txt")

rownames(counts) <- genes$V1
colnames(counts) <- cells$V1

counts <- counts[, grepl("10x_v2", colnames(counts))]

metadata <- read_tsv('data/meta.txt')
metadata <- filter(metadata, NAME %in% colnames(counts))

counts <- counts[, metadata$NAME]


sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata
)


# Sort gene names ---------------------------------------------------------

symbol <- sapply(strsplit(rownames(sce), "_"), `[`, 2)
is_dup <- duplicated(symbol)

sce <- sce[!is_dup,]
rownames(sce) <- symbol[!is_dup]


# QC ----------------------------------------------------------------------

is_mito <- grepl("^MT-", rownames(sce))

sce <- addPerCellQC(sce, subsets=list(mito=is_mito))

plotColData(sce, x="detected", y="subsets_mito_percent", colour_by = "CellType")

## Already QCd

sce <- computeSumFactors(sce)

sce <- logNormCounts(sce)

sce <- runPCA(sce, ncomponents = 50)
sce <- runUMAP(sce)


# Save --------------------------------------------------------------------

saveRDS(sce, "data/sce.rds")
