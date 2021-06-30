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

source("functions.R")

options(shiny.maxRequestSize=30*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("cytosel"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      fileInput("input_scrnaseq",
        "Input scRNA-seq",
        accept = c(".rds")
      ),
      selectInput(
          "coldata_column",
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

    # Show a plot of the generated distribution
    mainPanel(
      fluidRow(
        column(
          6,
          uiOutput("BL")
        ),
        column(
          6,
          plotOutput("all_plot"),
          plotOutput("top_plot")
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    umap_top <- reactiveVal()
    umap_all <- reactiveVal()
    column <- reactiveVal()
    
  output$all_plot <- renderPlot({
    req(umap_all)
      req(column)
      
      col <- column()
      
      if(is.null(umap_all())) { return(NULL) }

      ggplot(umap_all(), aes_string(x = "UMAP1", y = "UMAP2", color = col)) +
          geom_point() +
          labs(subtitle = "UMAP all genes")
  })
  
  output$top_plot <- renderPlot({
      req(umap_top)
      req(column)
      
      req(column)
      
      col <- column()
      
      if(is.null(umap_all())) { return(NULL) }
      
      ggplot(umap_top(), aes_string(x = "UMAP1", y = "UMAP2", color = col)) +
          geom_point() +
          labs(subtitle = "UMAP selected markers")
  })
  
  
  observeEvent(input$start_analysis, {
      req(input$coldata_column)
      req(input$panel_size)
      req(sce())
      
      column(input$coldata_column)
      
      withProgress(message = 'Processing data', value = 0, {
          
          incProgress(1, detail = "Finding markers")
        markers <- get_markers(sce(), column(), input$panel_size)
        
        output$BL <- renderUI({
            bucket_list(
            header = "Marker selection",
            orientation = "horizontal",
            group_name = "bucket_list_group",
            add_rank_list(
                text = "Selected markers",
                labels = markers$top_markers
            ),
            add_rank_list(
                text = "All genes/proteins",
                labels = markers$all_markers
            )
        )
        })
        
        incProgress(2, detail = "Computing UMAP")
        umap_all(get_umap(sce(), column()))
        umap_top(get_umap(sce()[markers$top_markers,], column()))
        
      })
      
      
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
  
}

# Run the application
shinyApp(ui = ui, server = server)
