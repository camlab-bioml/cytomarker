#' Make the expression and correlation heatmaps
#' 
#' @param sce A SingleCellExperiment object
#' @param markers A list storing three lists of markers: 
#' recommended_markers, scratch_markers, and top_markers
#' @param column A vector storing a column
#' @param display A vector storing selected heatmap display option
#' @param normalization A vector holding either 'Expression' or 'z-score'
#' @param pref_assay Assay loaded
#' 
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom stats cor
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly layout
create_heatmap <- function(sce, markers, column, display, normalization, pref_assay) {
  
  normalization <- match.arg(normalization, c("Expression", "z-score"))
  
  mat <- summarizeAssayByGroup(
    sce,
    id = colData(sce)[[column]],
    subset.row = markers$top_markers,
    statistics = 'mean',
    assay.type = pref_assay
  )
  
  mat <- (assay(mat, 'mean'))
  
  legend <- "Mean\nexpression"
  if(normalization == "z-score") {
    mat <- t(scale(t(mat)))
    legend <- "z-score\nexpression"
  }
  
  if(display == "Marker-marker correlation") {
    
    cc <- cor(t(mat))
    cor_map <- plot_ly(z=cc, 
                       type='heatmap',
                       x = rownames(cc),
                       y = colnames(cc)) %>% 
      layout(title='Correlation')
    
    return(cor_map)
  } else {
    
    
    
    expression_map <- plot_ly(z=t(mat), 
                              type='heatmap',
                              x = rownames(mat),
                              y = colnames(mat)) %>% 
      layout(title='Expression')
    
    return(expression_map)
  }
  
}

#' Return the legend of the cell type colours
#' 
#' @param palette A vector storing hex colours
#' 
#' @importFrom ggplot2 theme_bw geom_histogram scale_fill_manual guides
get_legend <- function(palette) {
  ## Draw the ggplot
  df <- tibble(r = stats::runif(length(palette)),
               `Cell type` = names(palette))
  
  # TODO: Maybe "Cell type" isn't the best name for this
  plt <- ggplot(df, aes(x = r, fill = `Cell type`)) + 
    geom_histogram() +
    scale_fill_manual(values = palette) +
    guides(fill = guide_legend(ncol = 8)) 
  
  legend <- cowplot::get_legend(plt)
  legend
}

#' Map gene symbols to antibody icons depending on
#' whether the gene is initially suggested or selected by the user
#' 
#' @param gene_id An element in the marker list
#' @param markers the list of markers
map_gene_name_to_antibody_icon <- function(gene_id, markers) {
  
  if(gene_id %in% markers$recommended_markers){
    return(icon("star"))
  }else{
    return(NULL)
  }
  
}