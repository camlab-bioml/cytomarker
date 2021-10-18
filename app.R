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
                  "Capture heterogeneity of:",
                  NULL,
                  multiple=TRUE),
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
      downloadButton('downloadData', 'Save\npanel'),
      width = 3
    ),
    mainPanel(tabsetPanel(
      tabPanel("Marker selection",
               div(style = "display:inline-block; horizontal-align:top; width:35%",
                   autocomplete_input("add_markers", "Add markers", options=c())),
               div(style = "display:inline-block; horizontal-align:top; width:35%",
                   actionButton("enter_marker", "Add marker")),
               fluidRow(column(12, uiOutput("BL")))),
      tabPanel("UMAP",
               fluidRow(
                 column(6, plotOutput("all_plot", height="600px")),
                 column(6, plotOutput("top_plot", height="600px"))
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
      tabPanel("Alternative markers",
               div(style = "display:inline-block; horizontal-align:top; width:35%",
                   autocomplete_input("input_gene", "Input gene", options=c())),
               div(style = "display:inline-block; horizontal-align:top; width:60%",
                   numericInput("number_correlations", "Number of alternative markers", 
                                value = 10, min = 1, 
                                width = "150px")),
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
    
    update_autocomplete_input(session, "input_gene",
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
      output$selected_assay <- renderText({paste("Assay type: ", pref_assay())})
    } else {
      showModal(assay_modal(failed = TRUE, input_assays()))
    }
  })

  output$all_plot <- renderPlot({
    req(umap_all)
    req(column)

    columns <- column()

    if (is.null(umap_all())) {
      return(NULL)
    }

    plts <- list()
    for(col in columns) {
      plts[[col]] <- ggplot(umap_all(), aes_string(x = "UMAP1", y = "UMAP2", color = col)) +
        geom_point() +
        labs(subtitle = "UMAP all genes")
    }
    cowplot::plot_grid(plotlist = plts, ncol=1)
  })
  
  output$top_plot <- renderPlot({
    req(umap_top)
    req(column)
    
    req(column)
    
    columns <- column()
    
    if (is.null(umap_top())) {
      return(NULL)
    }
    
    plts <- list()
    for(col in columns) {
      plts[[col]] <- ggplot(umap_top(), aes_string(x = "UMAP1", y = "UMAP2", color = col)) +
        geom_point() +
        labs(subtitle = "UMAP selected markers")
    }
    cowplot::plot_grid(plotlist = plts, ncol=1)
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
    mm <- metrics()
    if (is.null(mm)) {
      return(NULL)
    }
    
    columns <- names(mm)
    plots <- list()
    
    for(column in columns) {
      m <- mm[[column]]
    
      m$what <- fct_reorder(m$what, desc(m$score))
      m$what <- fct_relevel(m$what, "Overall")
      
      plots[[column]] <- ggplot(m, aes(y = fct_rev(what), x = score)) +
                    geom_boxplot(fill = 'grey90') +
                    labs(x = "Score",
                         y = "Source")
    }
    
    cowplot::plot_grid(plotlist = plots, ncol=1, labels = columns)
    
  })
  
  
  observeEvent(input$start_analysis, {
    req(input$coldata_column)
    req(input$panel_size)
    req(sce())
    
    ## Set initial markers:
    column(input$coldata_column)
    
    ## TODO: get markers for all columns, not just 1
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
      list(recommended_markers = input$bl_recommended,
           top_markers = input$bl_top,
           # all_markers = input$bl_all
           scratch_markers = input$bl_scratch)
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
  
  observeEvent(input$enter_marker, {
    
    if(!is.null(input$add_markers) && stringr::str_length(input$add_markers) > 1 && (input$add_markers %in% rownames(sce()))) {
      ## Need to update markers:
      new_marker <- input$add_markers
      
      cm <- current_markers()
      
      current_markers(
        list(recommended_markers = cm$recommended_markers,
             # all_markers = setdiff(cm$all_markers, new_marker),
             scratch_markers = setdiff(cm$scratch_markers, new_marker),
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
          text = "Recommended markers",
          labels = markers$recommended_markers,
          input_id = "bl_recommended",
          class = c("default-sortable", "cytocellbl")
        ),
        add_rank_list(
          text = "All genes/proteins",
          labels = markers$scratch_markers,
          # input_id = "bl_all",
          input_id = "bl_scratch",
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
      
      withProgress(message = 'Processing data', value = 0, {
        gene_to_replace(input$input_gene)
        
        incProgress(6, detail = "Computing alternatives")
      
        # Make this sampling dependent on the input sample argument
        x <- as.matrix(assay(sce(), pref_assay())[,sample(ncol(sce()), min(5000, ncol(sce())))])
        
        y <- x[gene_to_replace(), ]
        
        yo <- x[rownames(x) != gene_to_replace(),]
        
        correlations <- cor(t(yo), y)
        
        alternatives <- data.frame(Gene = rownames(yo), Correlation = correlations[,1])
        alternatives <- alternatives[!is.na(alternatives$Correlation),]
        alternatives <- alternatives[order(-(alternatives$Correlation)),]
        alternatives <- alternatives[1:input$number_correlations,]
        rownames(alternatives) <- NULL
        
        output$alternative_markers <- renderDT(alternatives, server = FALSE)
      })
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
      ## TODO: add multi column support to the following
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
