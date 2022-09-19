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
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
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
    
    cc <- round(cor(t(mat)), 4)
    cor_map <- plot_ly(z=cc, 
                       type='heatmap',
                       x = rownames(cc),
                       y = rownames(cc),
                       hovertemplate = "Gene(x): %{x}<br>Gene(y): %{y}<br>Correlation: %{z}<extra></extra>",
                       showticklabels = F) %>% 
      layout(title='Correlation', width = 650, height = 650)
    
    return(cor_map)
  } else {
    # Convert to dataframe to prevent issues when rownames are numbers
    df <- t(mat) |> 
      round(4) |> 
      as.data.frame() |> 
      rownames_to_column("y") |> 
      pivot_longer(-y, names_to = "gene", values_to = "expression")
    
    if(!all(is.na(as.numeric(df$y)))){
      # get row names and order (used to relevel factor)
      row_levels <- df$y |> 
        unique() |> 
        as.numeric() |> 
        sort()
      
      df <- mutate(df, y = factor(y, levels = row_levels))
    }
    
    expression_map <- plot_ly(df, x = ~gene, y = ~y, z = ~expression, 
                              type = "heatmap",
                              colorbar = list(title = as.character(normalization)),
                              hovertemplate = paste("Gene: %{x}<br>Cell Type: %{y}<br>",
                             as.character(normalization), ": %{z}<extra></extra>", sep = "")) |> 
      layout(title = as.character(normalization),
             xaxis = list(title = ''),
             yaxis = list(title = ''),
            height = 600, width = ifelse(25*length(unique(df$gene)) <= 450,
                                         450, ifelse(
                                           25*length(unique(df$gene)) <= 1000,
                                           25*length(unique(df$gene)), 1000
                                         )))
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
    guides(fill = guide_legend(ncol = 1)) 
  
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