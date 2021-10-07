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
library(readr)
library(dqshiny)
library(DT)

theme_set(theme_cowplot())

source("functions.R")

options(shiny.maxRequestSize = 1000 * 200 * 1024 ^ 2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("cytosel"),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "cytocell.css")
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("input_scrnaseq",
                "Input scRNA-seq",
                accept = c(".rds")),
      textOutput("selected_assay"),
      br(),
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
      checkboxInput("subsample_for_umap", "Subsample for UMAP", value = TRUE),
      actionButton("start_analysis", "Go"),
      actionButton("refresh_analysis", "Refresh"),
      hr(),
      downloadButton('downloadData', 'Download\nmarkers'),
      width = 3
    ),
    mainPanel(tabsetPanel(
      tabPanel("Marker selection",
               fluidRow(column(12, 
                               autocomplete_input("add_markers", "Add markers", options=c()),
                               )),
               fluidRow(column(12,
                              uiOutput("BL")))),
      tabPanel("UMAP",
               fluidRow(
                 column(6, plotOutput("all_plot")),
                 column(6, plotOutput("top_plot"))
               )),
      tabPanel(
        "Heatmap",
        selectInput(
          "heatmap_expression_norm",
          label = "Heatmap expression normalization",
          choices = c("Expression", "z-score")
        ),
        plotOutput("heatmap")
      ),
      tabPanel("Metrics",
               plotOutput("metric_plot")),
      tabPanel("Alternative Markers",
               textInput("input_gene", "Input gene"),
               actionButton("enter_gene", "Enter"),
               DTOutput("alternative_markers"))
    ))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  umap_top <- reactiveVal()
  umap_all <- reactiveVal()
  column <- reactiveVal()
  metrics <- reactiveVal()
  
  heatmap <- reactiveVal()
  current_markers <- reactiveVal()
  # heatmap_expression_norm <- reactiveVal()
  
  assay_modal <- function(failed = FALSE, assays) {
    modalDialog(
      # tags$script(HTML(enter_assay_selection)),
      selectInput("assay",
                  "Choose which assay to load",
                  assays),
      if (failed) {
        div(tags$b("Error", style = "color: red;"))
      },
      footer = tagList(
        actionButton("assay_select", "Select")
      )
    )
  }
  
  sce <- reactiveVal()
  input_assays <- reactiveVal()
  pref_assay <- reactiveVal("logcounts")
  
  observeEvent(input$input_scrnaseq, {
    sce(readRDS(input$input_scrnaseq$datapath))
    
    # Pop-up dialog box for assay selection
    input_assays(names(assays(sce())))
    showModal(assay_modal(assays = input_assays()))
    
    update_autocomplete_input(session, "add_markers",
                              options = rownames(sce()))
    
    updateSelectInput(
      session = session,
      inputId = "coldata_column",
      choices = colnames(colData(sce()))
    )
    
  })

  observeEvent(input$assay_select, {
    if (!is.null(input$assay)) {
      pref_assay(input$assay)
      removeModal()
    } else {
      showModal(assay_modal(failed = TRUE, input_assays()))
    }
  })
  
  output$selected_assay <- renderText({
    paste("Assay type: ", pref_assay())
  })

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
    if (is.null(m)) {
      return(NULL)
    }
    m$what <- fct_reorder(m$what, desc(m$score))
    m$what <- fct_relevel(m$what, "Overall")
    
    ggplot(m, aes(y = fct_rev(what), x = score)) +
      geom_boxplot(fill = 'grey90') +
      labs(x = "Score",
           y = "Cell type",
           subtitle = "Score for selected panel")
  })
  
  
  observeEvent(input$start_analysis, {
    req(input$coldata_column)
    req(input$panel_size)
    req(sce())
    
    ## Set initial markers:
    column(input$coldata_column)
    markers <- get_markers(sce(), column(), input$panel_size, pref_assay())
    current_markers(markers)
    
    update_analysis()
    
  })
  
  observeEvent(input$refresh_analysis, {
    req(input$coldata_column)
    req(input$panel_size)
    req(sce())
    
    ## Need to update markers:
    current_markers(
      list(top_markers = input$bl_top,
         all_markers = input$bl_all)
    )
    
    update_analysis()
    
  })
  
  observeEvent(input$heatmap_expression_norm, {
    ## What to do when heatmap selection is made
    req(sce())
    heatmap(
      create_heatmap(
        sce(),
        current_markers(),
        column(),
        input$heatmap_expression_norm,
        pref_assay()
      )
    )
  })
  
  observeEvent(input$add_markers, {
    new_marker <- input$add_markers
    if(!is.null(new_marker) && stringr::str_length(new_marker) > 1 && (new_marker %in% rownames(sce()))) {
      ## Need to update markers:
      print(new_marker)
      cm <- current_markers()
      
      current_markers(
        list(all_markers = setdiff(cm$all_markers, new_marker),
              top_markers = c(new_marker, setdiff(cm$top_markers, new_marker)))
      )
      
      update_BL(current_markers())
    }
  })
  
  update_BL <- function(markers) {
    output$BL <- renderUI({
      bucket_list(
        header = "Marker selection",
        orientation = "horizontal",
        group_name = "bucket_list_group",
        add_rank_list(
          text = "All genes/proteins",
          labels = markers$all_markers,
          input_id = "bl_all",
          class = c("default-sortable", "cytocellbl")
        ),
        add_rank_list(
          text = "Selected markers",
          labels = markers$top_markers,
          input_id = "bl_top",
          class = c("default-sortable", "cytocellbl")
        )
      )
    })
  }
  
  gene_to_replace <- reactiveVal()
  
  observeEvent(input$enter_gene, {
    req(input$input_gene)
    req(sce())
    
    if(input$input_gene %in% rownames(sce())) {
      gene_to_replace(input$input_gene)
      x <- assay(sce(), pref_assay())
      y <- x[gene_to_replace(),]
      yo <- x[rownames(x) != gene_to_replace(),]
      correlations <- cor(t(yo), y)
      
      alternatives <- data.frame(rownames(x), correlations)

      output$alternative_markers <- renderDT(alternatives[1:10,], server = FALSE) 
    }
    
  })
  
  update_analysis <- function() {
    column(input$coldata_column)
    
    withProgress(message = 'Processing data', value = 0, {
      markers <- current_markers()
      
      update_BL(markers)
      
      incProgress(2, detail = "Computing UMAP")
      
      if (input$subsample_for_umap) {
        nc <- ncol(sce())
        to_subsample <-
          sample(seq_len(nc), min(nc, 2000), replace = FALSE)
        umap_all(get_umap(sce()[, to_subsample], column()))
        umap_top(get_umap(sce()[markers$top_markers, to_subsample], column()))
      } else {
        umap_all(get_umap(sce(), column()))
        umap_top(get_umap(sce()[markers$top_markers,], column()))
      }
      
      incProgress(3, detail = "Drawing heatmap")
      heatmap(create_heatmap(sce(), markers, column(), input$heatmap_expression_norm))
      
      incProgress(4, detail = "Computing panel score")
      
      metrics(get_scores(sce(), column(), markers$top_markers, pref_assay()))
    })
  }
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('markers-', Sys.Date(), '.txt', sep='')
    },
    content = function(con) {
      selected_markers <- current_markers()$top_markers
      write_lines(selected_markers, con)
    }
  )
  
}




# Run the application
shinyApp(ui = ui, server = server)
