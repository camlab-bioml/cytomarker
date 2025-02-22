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
#' @importFrom stats cor
#' @importFrom viridis viridis
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
create_heatmap <- function(sce, markers, column, display, normalization, pref_assay) {
  
  library(heatmaply, quiet = T)
  
  normalization <- match.arg(normalization, c("Expression", "z-score"))
  
  mat <- scuttle::summarizeAssayByGroup(
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
    dims <- ifelse(length(rownames(cc)) < 70, 650, ifelse(length(rownames(cc)) < 95,
                                                          10.5*length(rownames(cc)), 1000))
    
    cor_map <- plotly::plot_ly(z=cc,
                       type='heatmap',
                       x = rownames(cc),
                       y = rownames(cc),
                       colorscale = list(c("blue"), c("red")),
                       hovertemplate = "Gene(x): %{x}<br>Gene(y): %{y}<br>Correlation: %{z}<extra></extra>",
                       showticklabels = T, width = dims,
                       height = dims) %>%
      plotly::layout(title='Correlation', yaxis = list(autotick = F, tickmode = "linear",
                                                       tickfont = list(size = ifelse(
                                                         length(rownames(cc)) > 75, 7.5, 10))),
             xaxis = list(autotick = F, tickmode = "linear", tickfont = list(size = ifelse(
               length(rownames(cc)) > 75, 7.5, 10))))
    # 
    # cor_map <- heatmaply::heatmaply(cc,
    #                      main = "Correlation",
    #                      plot_method = "plotly",
    #                      label_names = c("Gene", y = "Gene",
    #                                      "Correlation"),
    #                      key.title = "Correlation",
    #                      column_text_angle = 90) %>%
    #                       layout(width = 1.2*dims, height = dims)
    
    return(cor_map)
  } else {
    mat <- round(t(mat), 3)
    gene_num <- ncol(mat)
    
    expression_map <- heatmaply::heatmaply(mat,
                                main = as.character(normalization),
                                plot_method = "plotly",
                                label_names = c("Cell Type", y = "Gene", 
                                                as.character(normalization)),
                                key.title = as.character(normalization),
                                column_text_angle = 90,
                                height = 600, 
                                fontsize_col = ifelse(gene_num < 75, 10, ifelse(gene_num > 75 & 
                                                                                  gene_num <= 100, 7.5,
                                                                                5)),
                                width = ifelse(25*gene_num <= 450,
                                               450, ifelse(25*gene_num <= 1000,
                                                           25*gene_num, 1000)
                                )
    )
    
    return(expression_map)
  }
  
}

#' Return the legend of the cell type colours
#' 
# #' @param palette A vector storing hex colours
#' 
# #' @importFrom ggplot2 theme_bw geom_histogram scale_fill_manual guides
# get_legend <- function(palette) {
#   ## Draw the ggplot
#   df <- tibble(r = stats::runif(length(palette)),
#                `Cell type` = names(palette))
#   
#   # TODO: Maybe "Cell type" isn't the best name for this
#   plt <- ggplot(df, aes(x = r, fill = `Cell type`)) + 
#     geom_histogram() +
#     scale_fill_manual(values = palette) +
#     guides(fill = guide_legend(ncol = 1)) 
#   
#   legend <- cowplot::get_legend(plt)
#   
#   legend
# }

#' Map gene symbols to antibody icons depending on
#' whether the gene is initially suggested or selected by the user
#' 
#' @param gene_id An element in the marker list
#' @param markers the list of markers
map_gene_name_to_antibody_icon <- function(gene_id, markers) {
  
  # Option 1: if the gene is in the first run, label with a star
  # if(gene_id %in% markers$recommended_markers){
  if (gene_id %in% markers) {
    return(icon("star"))
  }else{
    return(NULL)
  }
  
}

#' Create a violin plot from a vector of genes, a cell type category in an SCE, and the assay desired
#' @param sce the SingleCellExperiment object
#' @param genes The vector of genes to visualize
#' @param cell_type The cell type category to resolve using the gene vector
#' @param assay The name of the assay to use. Common selections are logcounts or counts
make_violin_plot <- function(sce, genes, cell_type, assay) {
  expression_plot <- scater::plotExpression(
    sce, genes, x = cell_type,
    colour_by = cell_type,
    exprs_values = assay)
  
  return(as.data.frame(expression_plot["data"]) |> drop_na() |> 
           `colnames<-`(c("Gene", "Expression", "Cell Type", "Colour")) |>
           select(c(Gene, Expression, `Cell Type`, Colour)))
}
