utils::globalVariables(c("Cell type", "Symbol", "r", "score", "what"), "cytosel")

palette <- NULL


ggplot2::theme_set(cowplot::theme_cowplot())



## Global colour palette for cell types
palette <- NULL

options(shiny.maxRequestSize = 1000 * 200 * 1024 ^ 2)

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
#' @importFrom readr write_lines read_tsv
#' @importFrom dplyr desc
#' @importFrom bsplus use_bs_popover shinyInput_label_embed shiny_iconlink bs_embed_popover
#' @importFrom shinyjs useShinyjs hidden toggle
#' @importFrom grDevices dev.off pdf
#' @importFrom zip zip
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom randomcoloR distinctColorPalette
#' 
#' @param ... Additional arguments
cytosel <- function(...) {
  ui <- fluidPage(
    # Navigation prompt
    tags$head(
      tags$script(HTML("window.onbeforeunload = function() {return 'Your changes will be lost!';};"))
    ),
    
    # Use packages
    useShinyalert(),
    use_bs_popover(),
    useShinyjs(),
    
    # Title
    titlePanel("cytosel"),
    tags$head(
      includeCSS(system.file("www", "cytosel.css", package="cytosel"))
    ),
    
    # Side panel
    sidebarLayout(
      sidebarPanel(
        fileInput("input_scrnaseq", "Input scRNA-seq", accept = c(".rds")) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_embed_popover(content = get_tooltip('input_scrnaseq'),
                               placement = "right")
          ),
        textOutput("selected_assay"),
        br(),
        selectInput("coldata_column", "Capture heterogeneity of:", NULL, multiple=TRUE) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_embed_popover(content = get_tooltip('coldata_column'),
                      placement = "right", html = "true")
          ),
        numericInput("panel_size", "Targeted panel size:", 32, min = 1, max = 200, step = NA, width = NULL) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_embed_popover(content = get_tooltip('panel_size'),
                               placement = "right")
          ),
        radioButtons("marker_strategy", label = "Marker selection strategy",
                     choices = list("Cell type based"="fm", "Cell type free (geneBasis)" = "geneBasis"),
                     selected="fm"),
        checkboxInput("subsample_for_umap", "Subsample for UMAP", value = TRUE) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_embed_popover(content = get_tooltip('subsample_for_umap'),
                               placement = "right")
          ),
        actionButton("start_analysis", "Go"),
        actionButton("refresh_analysis", "Refresh"),
        hr(),
        downloadButton("downloadData", "Save\npanel"),
        width = 3
      ),
      
      # Tabs
      mainPanel(tabsetPanel(
        
        tabPanel("Marker selection",
                 icon = icon("list"),
                 br(),
                 fluidRow(column(5,
                          splitLayout(cellWidths = c(150, 150),
                                      div(style = "", autocomplete_input("add_markers", "Manual add markers", options=c(), width = "100%")),
                                      div(style = "margin-top:25px;", actionButton("enter_marker", "Add")))) ,
                          column(6,
                          splitLayout(cellWidths = c(200, 120, 140),
                                      div(style = "", fileInput("uploadMarkers", "Upload markers", width = "100%")),
                                      div(style = "margin-top:25px;", actionButton("add_to_selected", "Add uploaded", width = "100%")),
                                      div(style = "margin-top:25px;", actionButton("replace_selected", "Replace selected", width = "100%"))))),
                 hr(),
                 hidden(div(id = "marker_display",
                   tags$span(shiny_iconlink() %>%
                     bs_embed_popover(content = get_tooltip('marker_display'),
                                      placement = "top", html = TRUE)))),
                 plotOutput("legend", height='80px'),
                 hidden(div(id = "marker_visualization",
                   tags$span(shiny_iconlink() %>%
                     bs_embed_popover(content = get_tooltip('marker_visualization'),
                                      placement = "top", html = TRUE)))),
                 fluidRow(column(12, uiOutput("BL")))
          ),
        
        tabPanel("UMAP",
                 icon = icon("globe"),
                 br(),
                 helpText(get_tooltip('umap')),
                 fluidRow(column(6, plotOutput("all_plot", height="600px")),
                          column(6, plotOutput("top_plot", height="600px")))
          ),
        
        tabPanel("Heatmap",
                 icon = icon("table"),
                 br(),
                 splitLayout(cellWidths = c(320, 280),
                             div(selectInput("display_options", "Display expression or gene correlation", choices = c("Marker-marker correlation"), width = "86%") %>%
                                 shinyInput_label_embed(shiny_iconlink() %>%
                                                          bs_embed_popover(content = get_tooltip('heatmap_display_options'),
                                                                           placement = "right", html = "true"))),
                             hidden(div(id = "norm",
                                        selectInput("heatmap_expression_norm", "Heatmap expression normalization", choices = c("Expression", "z-score"), width = "89%") %>%
                                          shinyInput_label_embed(shiny_iconlink() %>%
                                                                   bs_embed_popover(content = get_tooltip('heatmap_expression_norm'),
                                                                                    placement = "right"))))),
                 plotOutput("heatmap"),
                 hr(),
                 tags$span(shiny_iconlink() %>%
                             bs_embed_popover(content = get_tooltip('gene_removal'),
                                              placement = "top")),
                 splitLayout(cellWidths = c(190, 200),
                             div(numericInput("n_genes", "Number of genes to remove", 
                                              value = 10, min = 1, width = "50%")),
                             div(style = "margin-top:25px;", actionButton("suggest_gene_removal", "View suggestions")))
          ),

        tabPanel("Metrics",
                 icon = icon("ruler"),
                 br(),
                 helpText(get_tooltip('metrics')),
                 tags$span(shiny_iconlink() %>%
                             bs_embed_popover(content = get_tooltip('metrics_explanation'),
                                              placement = "right", html = TRUE)),
                 textOutput("cells_per_category"),
                 plotOutput("metric_plot")
          ),

        tabPanel("Alternative markers",
                 icon = icon("exchange-alt"),
                 br(),
                 helpText(get_tooltip("alternative_markers")),
                 splitLayout(cellWidths = c(180, 240),
                             div(autocomplete_input("input_gene", "Input gene", options=c(), width = "100%")),
                             div(numericInput("number_correlations", "Number of alternative markers", value = 10, min = 1, width = "35%"))),
                 actionButton("enter_gene", "Enter"),
                 br(),
                 br(),
                 DTOutput("alternative_markers"),
                 hidden(div(id = "send", actionButton("send_markers", "Send markers to selection panel"))),
                 br(),
                 verbatimTextOutput("add_success")
          ),

        tabPanel("Antibody explorer",
                 icon = icon("wpexplorer"),
                 br(),
                 reactableOutput("antibody_table")
          )
        )
      )
    )
  )
  
  server <- function(input, output, session) {
    
    ### REACTIVE VARIABLES ###
    plots <- reactiveValues(all_plot = NULL, top_plot = NULL, metric_plot = NULL) # Save plots for download
    
    umap_top <- reactiveVal()
    umap_all <- reactiveVal()
    heatmap <- reactiveVal()
    metrics <- reactiveVal()
    
    sce <- reactiveVal()
    
    input_assays <- reactiveVal()
    pref_assay <- reactiveVal("logcounts")
    
    column <- reactiveVal()
    
    fms <- reactiveVal() # Where to store the findMarker outputs
    current_markers <- reactiveVal() # Currently selected markers
    num_selected <- reactiveVal()
    
    display <- reactiveVal() # Display options for heatmap
    
    replacements <- reactiveVal() # Where to store table of alternative markers
    
    suggestions <- reactiveVal() # Where to store suggested markers to remove
    
    cells_per_type <- reactiveVal() ## Number of cells per condition/type
    
    
    ### MODALS ###
    invalid_modal <- function() { # Input file is invalid
      shinyalert(
        title = "Error",
        text = paste("Input must be a Single Cell Experiment or Seurat object. Please upload a different file."),
        type = "error",
        showConfirmButton = TRUE,
        confirmButtonCol = "#337AB7"
      )
    }
    
    assay_modal <- function(failed = FALSE, assays) { # Assay selection
      modalDialog(
        selectInput("assay",
                    "Choose which assay to load",
                    assays),
        helpText("Recommended assay type is logcounts, as otherwise panel selection 
                 will be skewed towards high abundance transcripts rather than heterogeneously expressed transcripts."),
        if (failed) {
          div(tags$b("Error", style = "color: red;"))
        },
        footer = tagList(
          actionButton("assay_select", "Select")
        )
      )
    }
    
    unique_element_modal <- function(col) { # Column is invalid
      
      for(c in seq_len(length(col$colname))) {
        if(col$n[[c]] == 1) {
          shinyalert(
            title = "Error",
            text = paste("Only one level in column", col$colname[[c]], ". Please select another column."),
            type = "error",
            showConfirmButton = TRUE,
            confirmButtonCol = "#337AB7"
          )
        } else if(col$n[[c]] > 100) {
          shinyalert(
            title = "Error",
            text = paste("Column", col$colname[[c]], "has more than 100 unique elements. Please select another column."),
            type = "error",
            showConfirmButton = TRUE,
            confirmButtonCol = "#337AB7"
          )
        }
      }

    }
    
    warning_modal <- function(not_sce) { # Uploaded marker is not in SCE
      shinyalert(title = "Warning", 
                 text = paste(paste(not_sce, collapse = ", "), " are not in SCE."), 
                 type = "warning", 
                 showConfirmButton = TRUE,
                 confirmButtonCol = "#337AB7")
    }
    
    suggestion_modal <- function(failed = FALSE, suggestions) { # Marker removal suggestion
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
    
    dne_modal <- function(dne) { # Input marker is not in the dataset
      shinyalert(title = "Error", 
                 text = paste("Marker ", paste(dne), " does not exist in the dataset."), 
                 type = "error", 
                 showConfirmButton = TRUE,
                 confirmButtonCol = "#337AB7")
    }
    

    ### UPLOAD FILE ###
    observeEvent(input$input_scrnaseq, {
      
      obj <- readRDS(input$input_scrnaseq$datapath)
      
      if(isTruthy(methods::is(obj, 'SingleCellExperiment')) || isTruthy(methods::is(obj, 'Seurat'))) {
        
        sce(read_input_scrnaseq(obj))
        
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
        
      } else {
        invalid_modal()
      }
      
    })
    
    
    ### SELECT ASSAY TYPE ###
    observeEvent(input$assay_select, {
      if (!is.null(input$assay)) {
        pref_assay(input$assay)
        removeModal()
        output$selected_assay <- renderText({paste("Assay type: ", pref_assay())})
      } else {
        showModal(assay_modal(failed = TRUE, input_assays()))
      }
    })
    
    
    ### ANTIBODY EXPLORER ###
    output$antibody_table <- renderReactable({
      req(current_markers())
      
      markers <- current_markers()
      df_antibody <- dplyr::filter(antibody_info, Symbol %in% markers$top_markers)
      reactable(df_antibody,
                searchable = TRUE,
                filterable = TRUE,
                sortable = TRUE)
    })
    
    
    ### PLOTS ###
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
          labs(subtitle = "UMAP all genes") +
          scale_colour_manual(values=palette)
      }
      plots$all_plot <- cowplot::plot_grid(plotlist = plts, ncol=1)
      
      plots$all_plot
    },
    width=350,
    height=300)
    
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
          labs(subtitle = "UMAP selected markers") +
          scale_colour_manual(values=palette)
      }
      plots$top_plot <- cowplot::plot_grid(plotlist = plts, ncol=1)
      
      plots$top_plot
    },
      width=350,
      height=300
    )
    
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
      req(cells_per_type)
      mm <- metrics()
      if (is.null(mm)) {
        return(NULL)
      }
      
      columns <- names(mm)
      plts <- list()
      for(column in columns) {
        m <- mm[[column]]
        
        ## Add in number of cells per condition
        cpt <- cells_per_type()
        cpt['Overall'] <- sum(cpt)
        m$what <- plyr::mapvalues(as.character(m$what),
                                  from = names(cpt), 
                                  to = paste0(names(cpt), " (n = ", cpt, ")"))
        m$what <- as.factor(m$what)
        
        m$what <- fct_reorder(m$what, desc(m$score))
        m$what <- fct_relevel(m$what, paste0("Overall (n = ", cpt['Overall'], ")"))
        
        plts[[column]] <- ggplot(m, aes(y = fct_rev(what), x = score)) +
          geom_boxplot(fill = 'grey90') +
          labs(x = "Score",
               y = "Source")
      }
      plots$metric_plot <- cowplot::plot_grid(plotlist = plts, ncol=1, labels = columns)
      
      plots$metric_plot
    })
    
    
    ### ANALYSIS ###
    observeEvent(input$start_analysis, {

      withProgress(message = 'Initializing analysis', value = 0, {
        incProgress(detail = "Acquiring data")
        req(input$coldata_column)
        req(input$panel_size)
        req(sce())
        
        ## Set initial markers:
        scratch_markers_to_keep <- input$bl_scratch

        incProgress(detail = "Computing markers")
          
        columns <- good_col(sce(), input$coldata_column)
        column(columns$good)
        col <- columns$bad
        
        if(isTruthy(!is.null(column()))) {

          updateSelectInput(
            session = session,
            inputId = "display_options",
            choices = c("Marker-marker correlation", column())
          )
          
          if(!is.null(col)) {
            unique_element_modal(col)
          } 
          
          ## Get the markers first time          
          fms(
            compute_fm(sce(), column(), pref_assay())
          )
          
          markers <- get_markers(fms(), input$panel_size, input$marker_strategy, sce())
          markers$scratch_markers <- scratch_markers_to_keep
          
          # SMH
          current_markers(set_current_markers_safely(markers, fms()))
          
          num_selected(length(current_markers()$top_markers))
          cells_per_type(table(colData(sce())[[column()]]))
          
          update_analysis()
          
          
        } else {
          unique_element_modal(col)
        }
        
      })
      
    })
    
    observeEvent(input$refresh_analysis, {
      withProgress(message = 'Initializing analysis', value = 0, {
        req(column())
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
    
    
    ### MARKER SELECTION ###
    observeEvent(input$enter_marker, { # Manually add markers
      
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
        
        num_selected(length(current_markers()$top_markers))
        
        update_BL(current_markers(), num_selected())
        
      } else if(!(input$add_markers %in% rownames(sce()))) {
        dne_modal(dne = input$add_markers)
      }
    })
  
    observeEvent(input$add_to_selected, { # Add uploaded markers
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
    
    observeEvent(input$replace_selected, { # Replace selected markers by uploaded markers
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
    
    
    ### HEATMAP DISPLAY ###
    observeEvent(input$display_options, {
      display(input$display_options)
      
      toggle(id = "norm", condition = display() != "Marker-marker correlation")
      
      ## What to do when heatmap selection is made
      req(sce())
      
      columns <- column()
      
      observeEvent(input$heatmap_expression_norm, {
        
        if(display() == "Marker-marker correlation") { # Display heatmap for gene correlation
          
          for(col in columns) {
            heatmap(
              create_heatmap(sce(), current_markers(), col, display(), "Expression", pref_assay())
            )
          }
        } else { # Display heatmap for gene expression in a specific column
          
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
    
    
    ### REMOVE MARKERS ###
    observeEvent(input$suggest_gene_removal, { # Generate suggested markers to remove
      req(input$n_genes)
      
      expression <- as.matrix(assay(sce(), pref_assay())[current_markers()$top_markers,])
      cmat <- cor(t(expression))
      
      suggestions <- suggest_genes_to_remove(cmat, input$n_genes)
      
      showModal(suggestion_modal(suggestions = suggestions))
    })
    
    observeEvent(input$remove_suggested, { # Remove suggested markers
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
    
    
    ### ALTERNATIVE MARKERS ###
    observeEvent(input$enter_gene, { # Compute alternative markers
      req(input$number_correlations)
      req(sce())
      
      if(!is.null(input$input_gene) && stringr::str_length(input$input_gene) > 1 && (input$input_gene %in% rownames(sce()))) {
        
        withProgress(message = 'Processing data', value = 0, {
          incProgress(6, detail = "Computing alternatives")
          
          # Make this sampling dependent on the input sample argument
          replacements(
            compute_alternatives(input$input_gene, sce(), pref_assay(), input$number_correlations)
          )
          
          output$alternative_markers <- renderDT(replacements(), server = TRUE)
        })
        
      } else if(!(input$input_gene %in% rownames(sce()))) {
        dne_modal(dne = input$input_gene)
      }
    })
    
    observe({toggle(id = "send", condition = !is.null(input$alternative_markers_rows_selected) && !is.null(input$input_gene))})
      
    observeEvent(input$send_markers, { # Send alternative markers to selected markers panel
      n <- c(input$alternative_markers_rows_selected)
      
      replacements <- replacements()[n,]$Gene
      
      cm <- current_markers()
      markers <- list(recommended_markers = cm$recommended_markers,
                      scratch_markers = input$bl_scratch,
                      top_markers = unique(c(replacements, setdiff(cm$top_markers, input$bl_scratch))))
      
      # SMH
      current_markers(
        set_current_markers_safely(markers, fms())
      )
      
      num_selected(length(current_markers()$top_markers))

      update_BL(current_markers(), num_selected())
      
      output$add_success <- renderText({"Marker(s) added successfully."})
      
    })
    
    
    ### UPDATE SELECTED MARKERS ###
    update_BL <- function(markers, selected) {
      
      unique_cell_types <- sort(unique(markers$associated_cell_types))
      n_cell_types <- length(unique_cell_types)
      # palette <<- sample(cell_type_colors)[seq_len(n_cell_types)]
      set.seed(12345L)
      palette <<- distinctColorPalette(n_cell_types)
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
    
    
    ### UPDATE ANALYSIS ###
    update_analysis <- function() {
      
      withProgress(message = 'Processing data', value = 0, {
        markers <- current_markers()
        selected <- num_selected()
        
        update_BL(markers, selected)
        
        # Update UMAP
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
        
        # Update heatmap
        incProgress(detail = "Drawing heatmap")
        
        columns <- column()
        toggle(id = "norm", condition = display() != "Marker-marker correlation")
        
        if(display() == "Marker-marker correlation") {
          for(col in columns) {
            heatmap(
              create_heatmap(sce(), current_markers(), col, display(), "Expression", pref_assay())
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
        cells_per_type(table(colData(sce())[[column()]]))
        
        # Update metrics
        incProgress(detail = "Computing panel score")
        scores <- get_scores(sce(), column(), markers$top_markers, pref_assay())
        metrics(scores)
        
        
        # Show help text popover
        shinyjs::show(id = "marker_visualization")
        shinyjs::show(id = "marker_display")
        
      })
    }
    
    
    ### SAVE PANEL ###
    output$downloadData <- downloadHandler(
      filename = paste0("Cytosel-Panel-", Sys.Date(), ".zip"),
      content = function(fname) {
        download_data(fname, current_markers(), plots, heatmap(),
                      input_file = input$input_scrnaseq$datapath,
                      assay_used = pref_assay(),
                      het_source = column(),
                      panel_size = input$panel_size)
      },
      contentType = "application/zip"
    )
    
  }
  
  antibody_info <- read_tsv(system.file("inst", "abcam_antibodies_gene_symbol_associated.tsv", package="cytosel"))
  
  shinyApp(ui, server, ...)
}


