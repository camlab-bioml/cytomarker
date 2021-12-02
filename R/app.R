
ggplot2::theme_set(cowplot::theme_cowplot())



cell_type_colors <- c("#ED90A4", "#EC929A", "#E99590", "#E69787", "#E39A7D", "#DE9D74", "#D99F6B", "#D3A263", "#CDA55B",
"#C6A856", "#BEAB52", "#B6AE50", "#ADB150", "#A3B353", "#99B657", "#8EB85D", "#82BA65", "#76BC6D",
"#68BD76", "#59BF7F", "#48C089", "#33C192", "#14C19B", "#00C1A5", "#00C1AE", "#00C1B6", "#00C0BF",
"#00BFC7", "#00BDCE", "#19BCD5", "#39B9DB", "#50B7E0", "#63B4E4", "#75B0E8", "#85ADEA", "#94A9EC",
"#A2A5ED", "#AFA1EC", "#BA9EEB", "#C49AE8", "#CE97E5", "#D694E0", "#DC91DB", "#E290D5", "#E68ECE",
"#EA8EC7", "#EC8DBE", "#ED8EB6", "#EE8FAD", "#ED90A4")

## Global colour palette for cell types
palette <- NULL

# source("functions.R")

options(shiny.maxRequestSize = 1000 * 200 * 1024 ^ 2)

# Define UI for application that draws a histogram

#' Define main entrypoint of app
#' 
#' @export
#' 
#' @import shiny
#' @importFrom shinyalert useShinyalert shinyalert
#' @importFrom dqshiny autocomplete_input update_autocomplete_input
#' @importFrom DT DTOutput renderDT
#' @import SummarizedExperiment
#' @import ggplot2
#' @import forcats
#' @import sortable
#' @import reactable
#' @importFrom readr write_lines
#' @importFrom dplyr desc
#' @importFrom gridExtra grid.arrange
cytosel <- function(...) {
  ui <- fluidPage(
    # Navigation prompt
    tags$head(
      tags$script(HTML("window.onbeforeunload = function() {return 'Your changes will be lost!';};"))
    ),
    
    useShinyalert(),
    
    titlePanel("cytosel"),
    tags$head(
      includeCSS(system.file("www", "cytosel.css", package="cytosel"))
      # tags$link(rel = "stylesheet", type = "text/css", href = "cytosel.css")
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
                 icon = icon("list"),
                 # fluidRow(column(4, autocomplete_input("add_markers", "Manual add markers", options=c(),
                 #                                       width = "150px")),
                 #          column(6, div( style = "margin-top: 25px;", actionButton("enter_marker", "Add")))),
                 br(),
                 fluidRow(
                 column(6,
                 splitLayout(autocomplete_input("add_markers", "Manual add markers", options=c(),
                                                       width = "150px"),
                          div( style = "margin-top: 25px;", actionButton("enter_marker", "Add"))))
                 ),
                 hr(),
                 fluidRow(
                   column(12,
                     splitLayout(
                       div( style = "", fileInput("uploadMarkers", "Upload markers")),
                      div( style = "margin-top: 25px", actionButton("add_to_selected", "Add uploaded")),
                      div( style = "margin-top: 25px;", actionButton("replace_selected", "Replace selected")),
                      cellWidths = c(300, 120, 150)
                     )
                   )
                 ),
                 hr(),
                 plotOutput("legend", height='80px'),
                 fluidRow(column(12, uiOutput("BL")))
                 ),
        tabPanel("UMAP",
                 icon = icon("globe"),
                 fluidRow(column(6, plotOutput("all_plot", height="600px")),
                          column(6, plotOutput("top_plot", height="600px")))
                 ),
        tabPanel("Heatmap",
                 icon = icon("table"),
                 selectInput("display_options",
                             label = "Display heterogeneity expression or gene correlation",
                             choices = c()),
                 selectInput("heatmap_expression_norm",
                             label = "Heatmap expression normalization",
                             choices = c("Expression", "z-score")),
                 plotOutput("heatmap"),
                 hr(),
                 div(style = "display:inline-block; horizontal-align:top; width:20%",
                     numericInput("n_genes", "Number of genes to remove", 
                                  value = 10, min = 1,
                                  width = "190px")),
                 div(style = "display:inline-block; horizontal-align:top; width:75%",
                     actionButton("suggest_gene_removal", "View suggestions")),
                 # br(),
                 # textOutput("remove_gene"),
                 br()
                 ),
        tabPanel("Metrics",
                 icon = icon("ruler"),
                 plotOutput("metric_plot")
                 ),
        tabPanel("Alternative markers",
                 icon = icon("exchange-alt"),
                 div(style = "display:inline-block; horizontal-align:top; width:35%",
                     autocomplete_input("input_gene", "Input gene", options=c())),
                 div(style = "display:inline-block; horizontal-align:top; width:60%",
                     numericInput("number_correlations", "Number of alternative markers", 
                                  value = 10, min = 1, 
                                  width = "150px")),
                 actionButton("enter_gene", "Enter"),
                 DTOutput("alternative_markers")
                 ),
        tabPanel("Antibody explorer",
                 icon = icon("wpexplorer"),
                   reactableOutput("antibody_table")
                 )
      ))
    )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
    plots <- reactiveValues(all_plot = NULL, top_plot = NULL, metric_plot = NULL)
    
    umap_top <- reactiveVal()
    umap_all <- reactiveVal()
    column <- reactiveVal()
    metrics <- reactiveVal()
    
    heatmap <- reactiveVal()
    current_markers <- reactiveVal()
    num_selected <- reactiveVal()
    # heatmap_expression_norm <- reactiveVal()
    
    fms <- reactiveVal() # Where to store the findMarker outputs
    
    assay_modal <- function(failed = FALSE, assays) {
      modalDialog(
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
    
    unique_element_modal <- function(failed = FALSE, n_unique_elements, col) {
      modalDialog(
        textOutput("stop"),
        if (failed) {
          div(tags$b("Error", style = "color: red;"))
        },
        footer = tagList(
          modalButton("Cancel")
        )
      )
    }
    
    warning_modal <- function(not_sce) {
      shinyalert(title = "Warning", 
                 text = paste(paste(not_sce, collapse = ", "), " are not in SCE."), 
                 type = "warning", 
                 showConfirmButton = TRUE,
                 confirmButtonCol = "#337AB7")
    }
    
    suggestion_modal <- function(failed = FALSE, suggestions) {
      modalDialog(
        selectInput("markers_to_remove",
                    "Select markers to remove",
                    choices = current_markers()$top_markers,
                    selected = suggestions,
                    multiple = TRUE),
        if (failed) {
          div(tags$b("Error", style = "color: red;"))
        },
        footer = tagList(
          modalButton("Cancel"),
          actionButton("remove_suggested", "Move selected markers to scratch")
        )
      )
    }
    
    sce <- reactiveVal()
    input_assays <- reactiveVal()
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
    
    output$antibody_table <- renderReactable({
      req(current_markers())
      
      markers <- current_markers()
      df_antibody <- dplyr::filter(antibody_info, Symbol %in% markers$top_markers)
      reactable(df_antibody,
                searchable = TRUE,
                filter = TRUE,
                sortable = TRUE)
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
      plots$all_plot <- cowplot::plot_grid(plotlist = plts, ncol=1)
      
      plots$all_plot
    })
    
    output$top_plot <- renderPlot({
      req(umap_top)
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
      plots$top_plot <- cowplot::plot_grid(plotlist = plts, ncol=1)
      
      plots$top_plot
    })
    
    output$heatmap <- renderPlot({
      req(heatmap)
      req(column)
      
      if (is.null(heatmap())) {
        return(NULL)
      }

      ComplexHeatmap::draw(heatmap())

    })
    
    output$metric_plot <- renderPlot({
      req(metrics)
      mm <- metrics()
      if (is.null(mm)) {
        return(NULL)
      }
      
      columns <- names(mm)
      plts <- list()
      for(column in columns) {
        m <- mm[[column]]
      
        m$what <- fct_reorder(m$what, desc(m$score))
        m$what <- fct_relevel(m$what, "Overall")
        
        plts[[column]] <- ggplot(m, aes(y = fct_rev(what), x = score)) +
                      geom_boxplot(fill = 'grey90') +
                      labs(x = "Score",
                           y = "Source")
      }
      plots$metric_plot <- cowplot::plot_grid(plotlist = plts, ncol=1, labels = columns)
      
      plots$metric_plot
    })
    
    observeEvent(input$start_analysis, {
      
      withProgress(message = 'Initializing analysis', value = 0, {
        incProgress(detail = "Acquiring data")
        req(input$coldata_column)
        req(input$panel_size)
        req(sce())
        
        ## Set initial markers:
        column(input$coldata_column)
        
        columns <- column()
        
        updateSelectInput(
          session = session,
          inputId = "display_options",
          choices = c("Marker-marker correlation", columns)
        )
        
        scratch_markers_to_keep <- input$bl_scratch

        incProgress(detail = "Computing markers")
        
        ## Get the markers first time
        fms(
          compute_fm(sce(), columns, pref_assay())
        )
        
        markers <- get_markers(fms(), input$panel_size)
        markers$scratch_markers <- scratch_markers_to_keep
        
        # SMH
        current_markers(set_current_markers_safely(markers, fms()))
        
        num_selected(length(current_markers()$top_markers))
      
      })
      update_analysis()
    })
    
    observeEvent(input$refresh_analysis, {
      withProgress(message = 'Initializing analysis', value = 0, {
        req(input$coldata_column)
        req(input$panel_size)
        req(sce())
        
        incProgress(detail = "Recomputing markers")
        fms(
          compute_fm(sce(), column(), pref_assay())
        )
        # 
        markers <- list(recommended_markers = input$bl_recommended,
             scratch_markers = input$bl_scratch,
             top_markers = input$bl_top)
        
        # SMH
        current_markers(set_current_markers_safely(markers, fms()))
        
        num_selected(length(current_markers()$top_markers))
      })
      
      update_analysis()
    })
    
    display <- reactiveVal()
    
    observeEvent(input$display_options, {
      display(input$display_options)
      
      observeEvent(input$heatmap_expression_norm, {
        ## What to do when heatmap selection is made
        req(sce())
        
        columns <- column()
        
        if(display() == "Marker-marker correlation") {
          for(col in columns) {
            heatmap(
              create_heatmap(sce(), current_markers(), col, display(), input$heatmap_expression_norm, pref_assay())
            )
          }
        } else {
          for(col in columns) {
            if(col == display()) {
              heatmap(
                create_heatmap(sce(), current_markers(), display(), display(), input$heatmap_expression_norm, pref_assay())
              )
              break
            }
          }
        }
        
      })
      
    })
    
    observeEvent(input$enter_marker, {
      
      if(!is.null(input$add_markers) && stringr::str_length(input$add_markers) > 1 && (input$add_markers %in% rownames(sce()))) {
        ## Need to update markers:
        new_marker <- input$add_markers
        
        cm <- current_markers()
        markers <- list(recommended_markers = cm$recommended_markers,
                        scratch_markers = input$bl_scratch,
                        top_markers = unique(c(new_marker, setdiff(cm$top_markers, input$bl_scratch))))
        
        # SMH
        current_markers(
          set_current_markers_safely(markers, fms())
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
          
          ## Update markers
          markers <- list(recommended_markers = cm$recommended_markers,
                          scratch_markers = input$bl_scratch,
                          top_markers = unique(c(uploaded_markers[i], setdiff(cm$top_markers, input$bl_scratch))))
          
          # SMH
          current_markers(
            set_current_markers_safely(markers, fms())
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
          
          # SMH
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
    
    suggestions <- reactiveVal()
    
    observeEvent(input$suggest_gene_removal, {
      req(input$n_genes)
      
      expression <- as.matrix(assay(sce(), pref_assay())[current_markers()$top_markers,])
      cmat <- cor(t(expression))
      
      suggestions <- suggest_genes_to_remove(cmat, input$n_genes)
      
      showModal(suggestion_modal(suggestions = suggestions))
      # output$remove_gene <- renderText(paste("Suggested markers to remove: ", paste(suggestions, collapse = ', ')))
    })
    
    observeEvent(input$remove_suggested, {
      cm <- current_markers()
      
      if (!is.null(input$markers_to_remove)) {
        remove_marker <- c(input$markers_to_remove)
        num_selected(length(cm$top_markers))

        
        markers <- list(recommended_markers = cm$recommended_markers,
                   scratch_markers = unique(c(remove_marker, input$bl_scratch)),
                   top_markers = setdiff(cm$top_markers, remove_marker))       
        
        current_markers(
          set_current_markers_safely(markers, fms())
        )
        
        update_BL(current_markers(), num_selected())
        
        removeModal()
      } else {
        showModal(suggestion_modal(failed = TRUE, suggestions()))
      }

    })
    
    update_BL <- function(markers, selected) {
      
      unique_cell_types <- sort(unique(markers$associated_cell_types))
      n_cell_types <- length(unique_cell_types)
      set.seed(12345345L)
      palette <<- sample(cell_type_colors)[seq_len(n_cell_types)]
      names(palette) <<- unique_cell_types
      
      # markers$top_markers <- sapply(markers$top_markers, function(m) paste(icon("calendar"), m))
      
      labels_top <- lapply(markers$top_markers, 
                           function(x) div(map_gene_name_to_antibody_icon(x, antibody_info), x, style=paste('padding: 3px; background-color:', palette[ markers$associated_cell_types[x] ])))

      labels_recommended <- lapply(markers$recommended_markers, 
                           function(x) div(map_gene_name_to_antibody_icon(x, antibody_info), x, style=paste('padding: 3px; background-color:', palette[ markers$associated_cell_types[x] ])))

      labels_scratch <- lapply(markers$scratch_markers, 
                                   function(x) div(map_gene_name_to_antibody_icon(x, antibody_info), x, style=paste('padding: 3px; background-color:', palette[ markers$associated_cell_types[x] ])))
      
            
      output$legend <- renderPlot(cowplot::ggdraw(get_legend(palette)))
            
      output$BL <- renderUI({
        bucket_list(
          header = "Marker selection",
          orientation = "horizontal",
          group_name = "bucket_list_group",
          add_rank_list( 
            text = "Recommended markers",
            labels = labels_recommended,
            input_id = "bl_recommended",
            class = "cytocellbl",
            options = sortable_options(disabled = TRUE)
          ),
          add_rank_list(
            text = "Scratch space",
            labels = labels_scratch,
            input_id = "bl_scratch",
            class = c("default-sortable", "cytocellbl")
          ),
          add_rank_list(
            text = paste0("Selected markers (", selected, ")"),
            labels = labels_top,
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
        

        
        incProgress(detail = "Computing UMAP")
        
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
        
        incProgress(detail = "Drawing heatmap")
        columns <- column()
        if(display() == "Marker-marker correlation") {
          for(col in columns) {
            heatmap(
              create_heatmap(sce(), current_markers(), col, display(), input$heatmap_expression_norm, pref_assay())
            )
          }
        } else {
          for(col in columns) {
            if(col == display()) {
              heatmap(
                create_heatmap(sce(), current_markers(), display(), display(), input$heatmap_expression_norm, pref_assay())
              )
              break
            }
          }
        }
        
        incProgress(detail = "Computing panel score")
        scores <- get_scores(sce(), column(), markers$top_markers, pref_assay())
        metrics(scores)
        
      })
    }
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("Cytosel-", Sys.Date(), ".zip", sep = "")
      },
      content = function(fname) {
        
        path_marker_selection <- paste0("markers-", Sys.Date(), ".txt")
        path_umap <- paste0("UMAP-", Sys.Date(), ".pdf")
        path_heatmap <- paste0("heatmap-", Sys.Date(), ".pdf")
        path_metric <- paste0("metric-", Sys.Date(), ".pdf")
        
        selected_markers <- current_markers()$top_markers
        write_lines(selected_markers, path_marker_selection)
        
        pdf(path_umap, onefile = TRUE)
        grid.arrange(plots$all_plot, plots$top_plot)
        dev.off()
        
        pdf(path_heatmap, onefile = TRUE)
        draw(heatmap())
        dev.off()
        
        pdf(path_metric, onefile = TRUE)
        grid.arrange(plots$metric_plot)
        dev.off()
        
        fs <- c(path_marker_selection, path_umap, path_heatmap, path_metric)
        
        zip(zipfile = fname, files = fs) # zip function not working
      },
      contentType = "application/zip"
    )
    
  }
  
  antibody_info <- read_tsv(system.file("inst", "abcam_antibodies_gene_symbol_associated.tsv", package="cytosel"))
  
  shinyApp(ui, server, ...)
}


