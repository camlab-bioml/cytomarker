#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(sortable)
library(ggplot2)
library(scater)
library(forcats)
library(cowplot)

theme_set(theme_cowplot())

source("functions.R")

options(shiny.maxRequestSize = 1000 * 200 * 1024 ^ 2)

# Define UI for application that draws a histogram
ui <- fluidPage(titlePanel("cytosel"),
                tags$head(
                  tags$link(rel = "stylesheet", type = "text/css", href = "cytocell.css")
                ),
                sidebarLayout(
                  sidebarPanel(
                    fileInput("input_scrnaseq",
                              "Input scRNA-seq",
                              accept = c(".rds")),
                    selectInput("coldata_column",
                                "Cluster ID",
                                NULL),
                    numericInput(
                      "panel_size",
                      "Panel size",
                      32,
                      min = 1,
                      max = 200,
                      step = NA,
                      width = NULL
                    ),
                    actionButton("start_analysis", "Go"),
                    actionButton("refresh_analysis", "Refresh"),
                    width = 3
                  ),
                  mainPanel(tabsetPanel(
                    tabPanel("Marker selection", fluidRow(column(12,
                                                                 uiOutput(
                                                                   "BL"
                                                                 )))),
                    tabPanel("UMAP",
                             fluidRow(
                               column(6, plotOutput("all_plot")),
                               column(6, plotOutput("top_plot"))
                             )),
                    tabPanel("Heatmap", plotOutput("heatmap")),
                    tabPanel("Metrics",
                             plotOutput("metric_plot"))
                  ))
                ))

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  umap_top <- reactiveVal()
  umap_all <- reactiveVal()
  column <- reactiveVal()
  metrics <- reactiveVal()
  
  heatmap <- reactiveVal()
  
  output$all_plot <- renderPlot({
    req(umap_all)
    req(column)
    
    col <- column()
    
    if (is.null(umap_all())) {
      return(NULL)
    }
    
    ggplot(umap_all(), aes_string(x = "UMAP1", y = "UMAP2", color = col)) +
      geom_point() +
      labs(subtitle = "UMAP all genes")
  })
  
  output$top_plot <- renderPlot({
    req(umap_top)
    req(column)
    
    req(column)
    
    col <- column()
    
    if (is.null(umap_top())) {
      return(NULL)
    }
    
    ggplot(umap_top(), aes_string(x = "UMAP1", y = "UMAP2", color = col)) +
      geom_point() +
      labs(subtitle = "UMAP selected markers")
  })
  
  output$heatmap <- renderPlot({
    req(heatmap)
    req(column)

    col <- column()
    
    if (is.null(heatmap())) {
      return(NULL)
    }
    
    draw(heatmap())
  })
  
  output$metric_plot <- renderPlot({
    req(metrics)
    m <- metrics()
    if(is.null(m)) {
      return(NULL)
    }
    m$what <- fct_reorder(m$what, desc(m$score))
    m$what <- fct_relevel(m$what, "Overall")
    
    ggplot(m, aes(y = fct_rev(what), x = score)) +
      geom_boxplot(fill='grey90') +
      labs(x = "Score", y = "Cell type",
           subtitle = "Score for selected panel")
  })
  
  
  observeEvent(input$start_analysis, {
    req(input$coldata_column)
    req(input$panel_size)
    req(sce())
    
    update_analysis()
    
  })
  
  observeEvent(input$refresh_analysis, {
    req(input$coldata_column)
    req(input$panel_size)
    req(sce())
    
    update_analysis()
    
  })
  
  sce <- reactiveVal()
  
  observeEvent(input$input_scrnaseq, {
    print("Calling upload")
    print(input$input_scrnaseq$datapath)
    sce(readRDS(input$input_scrnaseq$datapath))
    
    
    updateSelectInput(
      session = session,
      inputId = "coldata_column",
      choices = colnames(colData(sce()))
    )
    
  })
  
  update_analysis <- function() {
    column(input$coldata_column)
    
    withProgress(message = 'Processing data', value = 0, {
      incProgress(1, detail = "Finding markers")
      markers <- get_markers(sce(), column(), input$panel_size)
      
      output$BL <- renderUI({
        bucket_list(
          header = "Marker selection",
          orientation = "horizontal",
          group_name = "bucket_list_group",
          add_rank_list(text = "All genes/proteins",
                        labels = markers$all_markers,
                        class=c("default-sortable", "cytocellbl")),
          add_rank_list(text = "Selected markers",
                        labels = markers$top_markers,
                        class=c("default-sortable", "cytocellbl"))
        )
      })
      
      incProgress(2, detail = "Computing UMAP")
      umap_all(get_umap(sce(), column()))
      umap_top(get_umap(sce()[markers$top_markers, ], column()))
      
      incProgress(3, detail = "Drawing heatmap")
      heatmap(create_heatmap(sce(), markers, column()))
      
      incProgress(4, detail = "Computing panel score")
      
      metrics(get_scores(sce(), column(), markers$top_markers))
      print(metrics())
      
    })
  }
  
}




# Run the application
shinyApp(ui = ui, server = server)
