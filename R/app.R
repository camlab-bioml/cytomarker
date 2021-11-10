
ggplot2::theme_set(cowplot::theme_cowplot())

# source("functions.R")

options(shiny.maxRequestSize = 1000 * 200 * 1024 ^ 2)

# Define UI for application that draws a histogram

#' Define main entrypoint of app
#' 
#' @import shiny
#' @importFrom shinyalert useShinyalert shinyalert
#' @importFrom dqshiny autocomplete_input update_autocomplete_input
#' @importFrom DT DTOutput renderDT
#' @import SummarizedExperiment
#' @import ggplot2
#' @import forcats
#' @import sortable
#' @importFrom readr write_lines
#' @importFrom ComplexHeatmap draw
#' @importFrom dplyr desc
cytosel <- function(...) {
  ui <- fluidPage(
    useShinyalert(),
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
        numericInput("panel_size",
                     "Targeted panel size:",
                     32, min = 1, max = 200,
                     step = NA, width = NULL),
        checkboxInput("subsample_for_umap", "Subsample for UMAP", value = TRUE),
        actionButton("start_analysis", "Go"),
        actionButton("refresh_analysis", "Refresh"),
        hr(),
        downloadButton("downloadData", "Save\npanel"),
        width = 3
      ),
      mainPanel(tabsetPanel(
        tabPanel("Marker selection",
                 div(style = "display:inline-block; horizontal-align:top; width:30%",
                     autocomplete_input("add_markers", "Add markers", options=c(),
                                         width = "150px")),
                 div(style = "display:inline-block; horizontal-align:top; width:65%",
                     fileInput("uploadMarkers", "Upload markers")),
                 br(),
                 div(style = "display:inline-block; horizontal-align:top; width:30%",
                     actionButton("enter_marker", "Add")),
                 div(style = "display:inline-block; horizontal-align:top; width:c(30%,30%)",
                     actionButton("add_to_selected", "Add"),
                     actionButton("replace_selected", "Replace all selected markers")),
                 hr(),
                 fluidRow(column(12, uiOutput("BL")))
                 ),
        tabPanel("UMAP",
                 fluidRow(column(6, plotOutput("all_plot", height="600px")),
                          column(6, plotOutput("top_plot", height="600px")))
                 ),
        tabPanel("Heatmap",
                 selectInput("heatmap_expression_norm",
                             label = "Heatmap expression normalization",
                             choices = c("Expression", "z-score")),
                 plotOutput("heatmap")
                 ),
        tabPanel("Metrics",
                 plotOutput("metric_plot")
                 ),
        tabPanel("Alternative markers",
                 div(style = "display:inline-block; horizontal-align:top; width:35%",
                     autocomplete_input("input_gene", "Input gene", options=c())),
                 div(style = "display:inline-block; horizontal-align:top; width:60%",
                     numericInput("number_correlations", "Number of alternative markers", 
                                  value = 10, min = 1, 
                                  width = "150px")),
                 actionButton("enter_gene", "Enter"),
                 DTOutput("alternative_markers")
                 )
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
    num_selected <- reactiveVal()
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
    
    warning_modal <- function(not_sce) {
      shinyalert(title = "Warning", 
                 text = paste(c(not_sce, " are not in SCE.")), 
                 type = "warning", 
                 showConfirmButton = TRUE,
                 confirmButtonCol = "#337AB7")
    }
    
    sce <- reactiveVal()
    pref_assay <- reactiveVal("logcounts")
    
    observeEvent(input$input_scrnaseq, {
      sce(read_input_scrnaseq(input$input_scrnaseq$datapath))
      
      input_assays <- c(names(assays(sce())))
      if("logcounts" %in% input_assays) {
        input_assays <- c("logcounts", input_assays[input_assays != "logcounts"])
      }
      
      showModal(assay_modal(assays = input_assays))
      
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
      
      columns <- column()
      
      scratch_markers_to_keep <- input$bl_scratch
      
      markers <- get_markers(sce(), columns, input$panel_size, pref_assay())
      markers$scratch_markers <- scratch_markers_to_keep
      
      current_markers(markers)
      
      num_selected(length(current_markers()$top_markers))
      
      update_analysis()
    })
    
    observeEvent(input$refresh_analysis, {
      req(input$coldata_column)
      req(input$panel_size)
      req(sce())
      
      ## Need to update markers:
      current_markers(
        list(recommended_markers = input$bl_recommended,
             scratch_markers = input$bl_scratch,
             top_markers = input$bl_top)
      )
      
      num_selected(length(current_markers()$top_markers))
      
      update_analysis()
    })
    
    observeEvent(input$heatmap_expression_norm, {
      ## What to do when heatmap selection is made
      req(sce())
      
      columns <- column()
      
      for(col in columns) {
        heatmap(
          create_heatmap(
            sce(),
            current_markers(),
            col,
            input$heatmap_expression_norm,
            pref_assay()
          )
        )
      }
    })
    
    observeEvent(input$enter_marker, {
      
      if(!is.null(input$add_markers) && stringr::str_length(input$add_markers) > 1 && (input$add_markers %in% rownames(sce()))) {
        ## Need to update markers:
        new_marker <- input$add_markers
        
        cm <- current_markers()
        
        current_markers(
          list(recommended_markers = cm$recommended_markers,
               scratch_markers = input$bl_scratch,
               top_markers = c(new_marker, setdiff(cm$top_markers, input$bl_scratch)))
        )
      }
      
      num_selected(length(current_markers()$top_markers))
        
      update_BL(current_markers(), num_selected())
    })
  
    observeEvent(input$add_to_selected, {
      req(input$uploadMarkers)
      
      uploaded_markers <- readLines(input$uploadMarkers$datapath)
      not_sce <- c()
      
      for(i in seq_len(length(uploaded_markers))) {
        if(!is.null(uploaded_markers[i]) && stringr::str_length(uploaded_markers[i]) > 1 && (uploaded_markers[i] %in% rownames(sce()))) {
          cm <- current_markers()
          current_markers(
            list(recommended_markers = cm$recommended_markers,
                 scratch_markers = input$bl_scratch,
                 top_markers = unique(c(uploaded_markers[i], setdiff(cm$top_markers, input$bl_scratch))))
          )
        } else if(!is.null(uploaded_markers[i]) && stringr::str_length(uploaded_markers[i]) > 1 && !(uploaded_markers[i] %in% rownames(sce()))) {
          not_sce <- c(not_sce, uploaded_markers[i])
        } 
      }
      
      if(length(not_sce) > 0) {
        warning_modal(not_sce)
      }
      
      num_selected(length(current_markers()$top_markers))
          
      update_BL(current_markers(), num_selected())
    })
    
    observeEvent(input$replace_selected, {
      req(input$uploadMarkers)
      
      uploaded_markers <- readLines(input$uploadMarkers$datapath)
      marker <- c()
      not_sce <- c()
      
      for(i in seq_len(length(uploaded_markers))) {
        if(!is.null(uploaded_markers[i]) && stringr::str_length(uploaded_markers[i]) > 1 && (uploaded_markers[i] %in% rownames(sce()))) {
          marker <- c(marker, uploaded_markers[i])
          
          cm <- current_markers()
          current_markers(
            list(recommended_markers = cm$recommended_markers,
                 scratch_markers = input$bl_scratch,
                 top_markers = marker)
          )
        } else if (!is.null(uploaded_markers[i]) && stringr::str_length(uploaded_markers[i]) > 1 && !(uploaded_markers[i] %in% rownames(sce()))) {
          not_sce <- c(not_sce, uploaded_markers[i])
        } 
      }
      
      if(length(not_sce) > 0) {
        warning_modal(not_sce)
      }
      
      num_selected(length(current_markers()$top_markers))
      
      update_BL(current_markers(), num_selected())
    })
    
    update_BL <- function(markers, selected) {
      output$BL <- renderUI({
        bucket_list(
          header = "Marker selection",
          orientation = "horizontal",
          group_name = "bucket_list_group",
          add_rank_list( 
            text = "Recommended markers",
            labels = markers$recommended_markers,
            input_id = "bl_recommended",
            class = "cytocellbl",
            options = sortable_options(disabled = TRUE)
          ),
          add_rank_list(
            text = "Scratch space",
            labels = markers$scratch_markers,
            input_id = "bl_scratch",
            class = c("default-sortable", "cytocellbl")
          ),
          add_rank_list(
            text = paste0("Selected markers (", selected, ")"),
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
        selected <- num_selected()
        
        update_BL(markers, selected)
        
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
        columns <- column()
        for(col in columns) {
          heatmap(create_heatmap(sce(), markers, col, input$heatmap_expression_norm, pref_assay()))
        } 
        
        incProgress(4, detail = "Computing panel score")
        scores <- get_scores(sce(), column(), markers$top_markers, pref_assay())
        metrics(scores)
        
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
  shinyApp(ui, server, ...)
}


