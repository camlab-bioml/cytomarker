
## Read the tooltips yaml on package load
tooltips <- yaml::read_yaml(
  # system.file("inst", "tooltips.yml", package="cytosel")
  file.path("inst", "tooltips.yml"), readLines.warn = F
)


#' Returns the tooltip as stored
#' in inst/tooltips.yml
get_tooltip <- function(tip_name) {
  tooltips[[tip_name]]
}

#' Get the tooltip info for the input SCE file
get_tooltip_q_input <- function() {
  return(paste('Upload SingleCellExperiment/Seurat data as an .rds file.',
  'Gene names should be in Gene Symbol format. If a dataset is too large (>1Gb),',
  'subsampling of cells is recommended. For some applications,',
  'filtering for protein coding-only genes may also increase the relevance ',
  'of the markers suggested by Cytosel.'))
}

#' Get the tooltip info for cell category input
get_tooltip_q_coldata_column <- function() {
 return(paste('Categories of interest must be created before using Cytosel',
              'Heterogeneity of multiple categories can be optimized at once if multiple categories are important.', 
              'Options include:',
                '<li>Cell type if there are known cell types that can be annotated.</li>',
                '<li>Cluster identity, if using clustering to use as a proxy for cell state / cell type.</li>',
                '<li>Condition, if aiming to separate experimental conditions.</li>',
                '<li>Timepoint, if aiming to separate distinct timepoints.</li>'))
}

