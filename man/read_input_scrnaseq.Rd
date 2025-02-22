% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{read_input_scrnaseq}
\alias{read_input_scrnaseq}
\title{Read the input single-cell RNA-seq dataset
from compressed RDS file. Accepts either:
(1) SingleCellExperiment
(2) Seurat object (converted to SingleCellExperiment)}
\usage{
read_input_scrnaseq(sce_path, filter_counts = T)
}
\arguments{
\item{sce_path}{Input uploaded path}

\item{filter_counts}{Boolean on whether to filter lowly expressed genes from the counts matrix(Default is True)}
}
\description{
Read the input single-cell RNA-seq dataset
from compressed RDS file. Accepts either:
(1) SingleCellExperiment
(2) Seurat object (converted to SingleCellExperiment)
}
