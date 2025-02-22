% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{get_umap}
\alias{get_umap}
\title{Compute the UMAP}
\usage{
get_umap(
  sce,
  columns,
  pref_assay,
  precomputed_vals = NULL,
  dim_col = NULL,
  only_top_markers = F,
  markers_to_use
)
}
\arguments{
\item{sce}{A SingleCellExperiment object}

\item{columns}{A vector storing columns}

\item{pref_assay}{Assay loaded}

\item{precomputed_vals}{Whether or not precomputed UMAP vals exist for the input SCE}

\item{dim_col}{If precomputed values exist, the name of the dimension holding the coordinates}

\item{only_top_markers}{Whether or not the UMAP computed will contain a subset of markers (Default is False)}

\item{markers_to_use}{The list of markers to compute the UMAP}
}
\description{
Compute the UMAP
}
