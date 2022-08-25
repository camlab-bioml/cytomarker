# utils::globalVariables(c("Cell type", "Symbol", "r", "score", "what"), "cytosel")

# palette <- NULL

utils::globalVariables(c("palette_to_use", "full_palette"), "cytosel")

ggplot2::theme_set(cowplot::theme_cowplot())

options(shiny.maxRequestSize = 1000 * 200 * 1024 ^ 2)

#' Define main entrypoint of app
#' 
#' @export
#' 
#' @import shiny shinytest
#' @importFrom shinyalert useShinyalert shinyalert
#' @importFrom dqshiny autocomplete_input update_autocomplete_input
#' @importFrom DT DTOutput renderDT
#' @importFrom tidyr drop_na
#' @import SummarizedExperiment
#' @import ggplot2
#' @import forcats
#' @import sortable
#' @import reactable
#' @importFrom readr write_lines read_tsv read_csv
#' @importFrom dplyr desc
#' @importFrom bsplus use_bs_popover shinyInput_label_embed shiny_iconlink bs_embed_popover
#' @importFrom shinyjs useShinyjs hidden toggle
#' @importFrom grDevices dev.off pdf
#' @importFrom zip zip
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom plotly plot_ly plotlyOutput renderPlotly layout
#' @importFrom parallelly availableCores
#' @importFrom BiocParallel MulticoreParam
#' @export
#' 
#' @param ... Additional arguments
cytosel <- function(...) {
  ## First up -- parse the antibodies
  # antibody_info <- read_tsv(system.file("inst", "abcam_antibodies_gene_symbol_associated.tsv", package="cytosel"))
  # antibody_info <- read_csv(system.file("inst", "Abcam_published_monos_with_gene.csv", package="cytosel"))
  
  antibody_info <- dplyr::rename(cytosel_data$antibody_info, Symbol = `Gene Name (Upper)`)
  antibody_info <- tidyr::drop_na(antibody_info)
  
  options(MulticoreParam=quote(MulticoreParam(workers=availableCores())))
  
  
  applications_parsed <- get_antibody_applications(antibody_info, 
                                                   'Symbol', 'Listed Applications')
  
  grch38 <- cytosel_data$grch38
  
  full_palette <- create_global_colour_palette()
  
  ui <- fluidPage(
    # Navigation prompt
    tags$head(
      tags$script(HTML("window.onbeforeunload = function() {return 'Your changes will be lost!';};"))
    ),
    tags$head(tags$style(".modal-dialog{ width:750px}")),
    
    # Use packages
    useShinyalert(force = TRUE),
    use_bs_popover(),
    useShinyjs(),
    
    # Title
    titlePanel("cytosel"),
    tags$head(
      # includeCSS(system.file("www", "cytosel.css", package="cytosel"))
      includeCSS(system.file(file.path("www", "cytosel.css"),
                             package = "cytosel"))
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
        selectInput("coldata_column", "Capture heterogeneity of:", NULL, multiple=FALSE) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_embed_popover(content = get_tooltip('coldata_column'),
                      placement = "right", html = "true")
          ),
        actionButton("show_cat_table", "Summary/Subset selection",
                     align = "center",
                     width = '100%', style='font-size:85%'),
        numericInput("panel_size", "Targeted panel size:", 32, min = 1, max = 200, step = NA, width = NULL) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_embed_popover(content = get_tooltip('panel_size'),
                               placement = "right")
          ),
        numericInput("min_category_count", "Minimum cell category cutoff:", 2, min = 2, max = 100, step = 0.5, width = NULL) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_embed_popover(content = get_tooltip('cat_cutoff'),
                               placement = "right")
          ),
        radioButtons("marker_strategy", label = "Marker selection strategy",
                     choices = list("Cell type based"="fm", "Cell type free (geneBasis)" = "geneBasis"),
                     selected="fm"),
        checkboxInput("subsample_sce", "Subsample cells", value = TRUE) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_embed_popover(content = get_tooltip('subsample_sce'),
                               placement = "right")
          ),
        br(),
        selectInput("select_aa", "Antibody applications:", applications_parsed$unique_applications, multiple=TRUE) %>%
          shinyInput_label_embed(
            shiny_iconlink() %>%
              bs_embed_popover(content = "Placeholder",
                               placement = "right", html = "true")
          ),
        
        actionButton("start_analysis", "Go", icon=icon("play")),
        # actionButton("refresh_analysis", "Refresh"),
        hr(),
        downloadButton("downloadData", "Save\npanel"),
        width = 3
      ),
      
      # Tabs
      mainPanel(tabsetPanel(id = "analysis_panels",
        tabPanel("Marker selection",
                 icon = icon("list"),
                 br(),
                 fluidRow(column(5,
                          splitLayout(cellWidths = c(150, 150),
                                      div(style = "", autocomplete_input("add_markers", "Manual add markers", options=c(), width = "100%",
                                                                         max_options = 6, contains = T)),
                                      div(style = "margin-top:25px;", actionButton("enter_marker", "Add")))) ,
                          column(6,
                          splitLayout(cellWidths = c(200, 120, 140),
                                      div(style = "", fileInput("uploadMarkers", "Upload markers", width = "100%")),
                                      div(style = "margin-top:25px;", actionButton("add_to_selected", "Add uploaded", width = "100%")),
                                      div(style = "margin-top:25px;", actionButton("replace_selected", "Replace selected", width = "100%"))))),
                 fluidRow(column(12, splitLayout(cellWidths = c(150, 150),
                                                selectInput('cell_type_markers', "Suggest markers for cell type:", choices=NULL),
                                            div(style = "margin-top:25px;", actionButton('add_cell_type_markers', "Suggest"),
                                                )))
                          ),
                 hr(),
                 hidden(div(id = "marker_display",
                   tags$span(shiny_iconlink() %>%
                     bs_embed_popover(content = get_tooltip('marker_display'),
                                      placement = "top", html = TRUE)))),
                 plotOutput("legend", height='80px'),
                 fluidRow(column(1, actionButton("refresh_marker_counts", "Refresh Marker Space counts",
                                                 style='font-size:85%'),
                                 align = "center"
                                 , style = "margin-bottom: 10px;"
                                 , style = "margin-top: 10px;")),
                 hidden(div(id = "marker_visualization",
                   tags$span(shiny_iconlink() %>%
                     bs_embed_popover(content = get_tooltip('marker_visualization'),
                                      placement = "top", html = TRUE)))),
                 fluidRow(column(12, uiOutput("BL"),
                                 ))
          ),
        
        tabPanel("UMAP",
                 icon = icon("globe"),
                 br(),
                 helpText(get_tooltip('umap')),
                 fluidRow(column(6, plotlyOutput("all_plot", width="500px", height="350px")),
                          column(6, plotlyOutput("top_plot", width="500px", height="350px")))
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
                 plotlyOutput("heatmap", height="600px"),
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
                 plotlyOutput("metric_plot")
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
                 br()
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
    allowed_genes <- reactiveVal() # Set of "allowed" genes (those with anotibodies, not ribo or mito)
    current_markers <- reactiveVal() # Currently selected markers
    
    display <- reactiveVal() # Display options for heatmap
    
    replacements <- reactiveVal() # Where to store table of alternative markers
    
    suggestions <- reactiveVal() # Where to store suggested markers to remove
    
    cells_per_type <- reactiveVal() ## Number of cells per condition/type
    
    summary_cat_tally <- reactiveVal()
    
    cytosel_palette <- reactiveVal()
    
    proceed_with_analysis <- reactiveVal(TRUE)
    
    specific_cell_types_selected <- reactiveVal()
    
    cell_types_high_enough <- reactiveVal()
    
    cell_types_to_keep <- reactiveVal()
    
    ## Current data frame of selected cell type markers and the reactable
    current_cell_type_marker_fm <- NULL
    selected_cell_type_markers <- reactive(getReactableState("cell_type_marker_reactable", "selected"))
    
    any_cells_present <- reactiveVal(TRUE)
    
    num_markers_in_selected <- reactiveVal()
    num_markers_in_scratch <- reactiveVal()
    
    cell_types_excluded <- reactiveVal()
    
    marker_suggestions <- reactiveVal()
    alternative_marks <- reactiveVal()

    
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
    
    assay_modal <- function(assays, failed = FALSE) { # Assay selection
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
                    selected = suggestions(),
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
                 text = paste("Marker ", paste(dne), " does not exist or has no corresponding antibody."), 
                 type = "error", 
                 showConfirmButton = TRUE,
                 confirmButtonCol = "#337AB7")
    }
    
    cell_cat_modal <- function() {
      modalDialog(
        DT::dataTableOutput("cell_cat_table"),
        selectInput("user_selected_cells", "Create a custom subset for analysis", NULL, multiple=TRUE),
        title = "Frequency Count for selected heterogeneity category",
        helpText("Certain cell types can be manually ignored by the user during analysis in the dialog box above. If this cell is left empty, then by default 
                 all cell types with a Freq of 2 or greater will be retained for analysis."),
        size = "xl",
        easyClose = TRUE,        
        footer = tagList(
          actionButton("add_selected_to_analysis", "Set subset for analysis."),
          modalButton("Exit.")
        ))
    }
    
    threshold_too_low_modal <- function() { # Input marker is not in the dataset
      shinyalert(title = "Error", 
                 text = "Cytosel requires the minimum cell type category to be set to 2 or greater for statistical inference. Please adjust the value of minimum cell category cutoff to at least 2.", 
                 type = "error", 
                 showConfirmButton = TRUE,
                 confirmButtonCol = "#337AB7")
    }

    cell_type_ignored_modal <- function() {
      
      shinyalert(title = "Warning",
                 text = HTML(paste("Cell types with abundance below the set threshold of ",
                                   input$min_category_count,
                                      ":", '<br/>',
                                   "<b>", toString(cell_types_excluded()), "</b>", '<br/>',
                                      "were identified and not removed by the user.",
                                      "They are to be ignored during analysis. The threshold can be changed using Minimum cell category cutoff.")),
                 type = "warning",
                 showConfirmButton = TRUE,
                 confirmButtonCol = "#337AB7",
                 html = TRUE)
    }
    
    no_cells_left_modal <- function() { # Uploaded marker is not in SCE
      shinyalert(title = "Error", 
                 text = "No cells remain after filtering/subsetting. Please review the input parameters.", 
                 type = "error", 
                 showConfirmButton = TRUE,
                 confirmButtonCol = "#337AB7")
    }
    
    
    ### UPLOAD FILE ###
    observeEvent(input$input_scrnaseq, {
      #if(isTruthy(methods::is(obj, 'SingleCellExperiment')) || isTruthy(methods::is(obj, 'Seurat'))) {
      updateNumericInput(session, "min_category_count", value = 2)
      input_sce <- read_input_scrnaseq(input$input_scrnaseq$datapath)
      input_sce <- detect_assay_and_create_logcounts(input_sce)
      input_sce <- parse_gene_names(input_sce, grch38)
      sce(input_sce)
      input_assays <- c(names(assays(sce())))

      # If there is more than 1 assay user to select appropriate assay
      if(length(input_assays) > 1){
        if("logcounts" %in% input_assays) {
          input_assays <- c("logcounts", input_assays[input_assays != "logcounts"])
        }

        showModal(assay_modal(assays = input_assays))
      }else{
        throw_error_or_warning(message = paste("Only one assay provided, thus using",
                                               input_assays),
                               duration = 5,
                               notificationType = 'message')
      }

      updateSelectInput(
        session = session,
        inputId = "coldata_column",
        choices = colnames(colData(sce()))
      )

    })
    
    observeEvent(input$coldata_column, {
      req(input$input_scrnaseq)
      specific_cell_types_selected(NULL)
      })
    
    observeEvent(input$add_selected_to_analysis, {
      req(input$input_scrnaseq)
      req(input$coldata_column)
      specific_cell_types_selected(input$user_selected_cells)
      
      if (length(specific_cell_types_selected()) > 0) {
        showNotification("Setting subset to select cell types",
                         duration = 2)
      } else {
        showNotification("Empty subset selection, defaulting to using all cell types in analysis.",
                         duration = 2)
      }
      
    })
    
    observeEvent(input$show_cat_table, {
    req(input$input_scrnaseq)
    req(input$assay_select)
    req(input$coldata_column)
    req(input$min_category_count)
    req(sce())
    
    cell_category_table <- create_table_of_hetero_cat(sce(),
                                                        input$coldata_column)
    
    if (nrow(cell_category_table) > 100) {
      shinyalert(
        title = "Error",
        text = paste("Column", input$coldata_column, "has more than 100 unique elements. Please select another column."),
        type = "error",
        showConfirmButton = TRUE,
        confirmButtonCol = "#337AB7"
      )
    } else {
      summary_cat_tally(cell_category_table)
      
      output$cell_cat_table <- DT::renderDataTable(
        summary_cat_tally()
      )
      
      if (!isTruthy(specific_cell_types_selected())) {
        specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      }
      
      updateSelectInput(
        session,
        "user_selected_cells",
        choices = unique(sce()[[input$coldata_column]]),
        selected = specific_cell_types_selected()
      )
      
      showModal(cell_cat_modal())
    }
    })
    
    ### SELECT ASSAY TYPE ###
    observeEvent(input$assay_select, {
      if (!is.null(input$assay)) {
        pref_assay(input$assay)
        removeModal()
        output$selected_assay <- renderText({paste("Selected assay: ", pref_assay())})
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
    output$all_plot <- renderPlotly({
      req(umap_all)
      req(column)
      
      columns <- column()
      
      if (is.null(umap_all())) {
        return(NULL)
      }
      
      # plts <- list()
      # for(col in columns) {
        # plts[[col]] <- ggplot(umap_all(), aes_string(x = "UMAP1", y = "UMAP2", color = col)) +
        #   geom_point() +
        #   labs(subtitle = "UMAP all genes") +
        #   scale_colour_manual(values=palette)
      # }
      # plots$all_plot <- cowplot::plot_grid(plotlist = plts, ncol=1)
      # plots$all_plot <- 
      
      # plots$all_plot

      plot_ly(umap_all(), x=~UMAP1, y=~UMAP2, color=~get(columns[1]), text=~get(columns[1]), 
              type='scatter', hoverinfo="text", colors=cytosel_palette()) %>% 
        layout(title = "UMAP all genes")
    })
    # },
    # width=350,
    # height=300)
    
    output$top_plot <- renderPlotly({
      req(umap_top)
      req(column)
      
      columns <- column()
      
      if (is.null(umap_top())) {
        return(NULL)
      }
      
      # plts <- list()
      # for(col in columns) {
      #   plts[[col]] <- ggplot(umap_top(), aes_string(x = "UMAP1", y = "UMAP2", color = col)) +
      #     geom_point() +
      #     labs(subtitle = "UMAP selected markers") +
      #     scale_colour_manual(values=palette)
      # }
      # plots$top_plot <- cowplot::plot_grid(plotlist = plts, ncol=1)
      # 
      # plots$top_plot

      plot_ly(umap_top(), x=~UMAP1, y=~UMAP2, color=~get(columns[1]), text=~get(columns[1]), 
              type='scatter', hoverinfo="text", colors=cytosel_palette()) %>% 
        layout(title = "UMAP selected markers")
    })
    
    output$heatmap <- renderPlotly({
      req(heatmap)
      req(column)
      
      if (is.null(heatmap())) {
        return(NULL)
      }
      
      return(heatmap())
    })
    
    output$metric_plot <- renderPlotly({
      req(metrics)
      req(cells_per_type)
      mm <- metrics()
      if (is.null(mm)) {
        return(NULL)
      }
      columns <- names(mm)
      plts <- list()
      column <- columns[1]
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
      m$what <- fct_rev(m$what)
      
      plot_ly(m, x = ~score, y = ~what, type='box', hoverinfo = 'none') %>% 
        layout(xaxis = list(title="Score"),
               yaxis = list(title="Source"))
      
    })
    
    ### ANALYSIS ###
    observeEvent(input$start_analysis, {
      req(input$coldata_column)
      req(input$panel_size)
      req(input$min_category_count)
      req(sce())
      
      if (input$min_category_count < 2) {
        proceed_with_analysis(FALSE)
        threshold_too_low_modal()
      } else {
        proceed_with_analysis(TRUE)
      }
      
      # if the user never opened up the tally table, automatically set all cell types
      # for analysis
      if (!isTruthy(specific_cell_types_selected())) {
        specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      }
      
      withProgress(message = 'Initializing analysis', value = 0, {
        incProgress(detail = "Checking data")
        
        cell_cat_summary <- create_table_of_hetero_cat(sce(),
                                                       input$coldata_column)
        
        cell_types_high_enough(remove_cell_types_by_min_counts(cell_cat_summary,
                                                               sce(), 
                                                               input$coldata_column, 
                                                               input$min_category_count))
        
        cell_types_to_keep(intersect(specific_cell_types_selected(),
                                     cell_types_high_enough()))
        
        sce(create_sce_column_for_analysis(sce(), cell_types_to_keep(), 
                                           input$coldata_column))
        
        if (dim(sce()[,sce()$keep_for_analysis == "Yes"])[2] == 0) {
          no_cells_left_modal()
          any_cells_present(FALSE)
        } else {
          any_cells_present(TRUE)
        }
        
        different_cells <- setdiff(specific_cell_types_selected(),
                                   cell_types_high_enough())
        
        if (length(different_cells) > 0 & isTruthy(any_cells_present())) {
          output$ignored_cell_types <- renderDataTable(data.frame(
            `Cell Type Ignored` = different_cells))
          cell_types_excluded(different_cells[!is.null(different_cells)])
          cell_type_ignored_modal()
        }
        
      })
      
      withProgress(message = 'Conducting analysis', value = 0, {
        incProgress(detail = "Acquiring data")
        req(proceed_with_analysis())
        req(any_cells_present())
          
          ## Set initial markers:
          scratch_markers_to_keep <- input$bl_scratch
          
          incProgress(detail = "Computing markers")
          
          columns <- good_col(sce()[,sce()$keep_for_analysis == "Yes"], input$coldata_column)
          column(columns$good)
          col <- columns$bad
          
          
          if (input$subsample_sce) {
            incProgress(detail = "Subsampling data")
            
            avail_cols <- which(sce()$keep_for_analysis == "Yes")
            to_subsample <- sample(avail_cols, min(length(avail_cols), 2000),
                                   replace = FALSE)
            
            sce(create_keep_vector_during_subsetting(sce(), to_subsample))
            
          } 
          
          if(isTruthy(!is.null(column()))) {
            
            updateSelectInput(
              session = session,
              inputId = "display_options",
              choices = c("Marker-marker correlation", column())
            )
            
            if(!is.null(col)) {
              unique_element_modal(col)
            } 
            
            ## Compute set of allowed genes
            allowed_genes(
              get_allowed_genes(input$select_aa, applications_parsed, 
                                sce()[,sce()$keep_for_analysis == "Yes"])
            )
            
            
            ## Change selected column to character to avoid factor levels without data
            # sce <- sce()
            sce(convert_column_to_character_or_factor(sce(), column()))
            # sce(sce)
            
            ## Get the markers first time 
            
            fms(
              compute_fm(sce()[,sce()$keep_for_analysis == "Yes"], 
                         column(), 
                         pref_assay(),
                         allowed_genes()
              )
            )
            
            updateSelectInput(
              session = session,
              inputId = "cell_type_markers",
              choices = names(fms()[[1]])
            )
            
            if(!is.null(input$bl_top)) {
              ## We get here if input$bl_top exists, ie if this
              ## is an analysis refresh
              ## in this case we set the markers to their existing values
              markers <- list(recommended_markers = input$bl_recommended,
                              scratch_markers = input$bl_scratch,
                              top_markers = input$bl_top)
            } else {
              ## We compute the set of markers for the first time
              markers <- get_markers(fms(), 
                                     # Adding 10 to make sure panel size is approximate 
                                     # since a) the same marker is selected multiple times and
                                     # b) excess markers are removed
                                     input$panel_size, 
                                     input$marker_strategy, 
                                     sce()[,sce()$keep_for_analysis == "Yes"],
                                     allowed_genes())
              
              if(length(markers$recommended_markers) < input$panel_size){
                showNotification("The cell types of the uploaded dataset show expression redundancy.\n
                               This results in fewer genes being shown than requested.",
                                 type = 'message',
                                 duration = NULL)
              }
              ## Forgotten what this is for
              markers$scratch_markers <- scratch_markers_to_keep
            }
            
            # SMH
            current_markers(set_current_markers_safely(markers, fms()))
          
            num_markers_in_selected(length(current_markers()$top_markers))
            num_markers_in_scratch(length(current_markers()$scratch_markers))
            cells_per_type(table(colData(
              sce()[,sce()$keep_for_analysis == "Yes"])[[column()]]))
            
            update_analysis()
            
            
          } else {
            unique_element_modal(col)
          }
          
        
      })
      
    })
    
    # # re-count the number of markers in each space when the sortable js is changed
    # observeEvent(input$bl_scratch, {
    #   markers <- list(recommended_markers = current_markers()$recommended_markers,
    #                   scratch_markers = input$bl_scratch,
    #                   top_markers = current_markers()$top_markers)
    #   current_markers(set_current_markers_safely(markers, fms()))
    #   num_markers_in_selected(length(current_markers()$top_markers))
    #   num_markers_in_scratch(length(current_markers()$scratch_markers))
    #   update_BL(current_markers(), num_markers_in_selected(),
    #             num_markers_in_scratch(),
    #             names(fms()[[1]]))
    # })
    
    # re-count the number of markers in each space when the sortable js is changed
    observeEvent(input$refresh_marker_counts, {
      req(sce())
      req(current_markers())
      
      current_markers(list(recommended_markers = current_markers()$recommended_markers,
                      scratch_markers = input$bl_scratch,
                      top_markers = input$bl_top,
                      associated_cell_types = current_markers()$associated_cell_types))
      
      # current_markers(set_current_markers_safely(markers, fms()))
      num_markers_in_selected(length(current_markers()$top_markers))
      num_markers_in_scratch(length(current_markers()$scratch_markers))
      unique_cell_types <- unique(unlist(current_markers()$associated_cell_types))
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                unique_cell_types)
    })
    
    
    # observeEvent(input$refresh_analysis, {
    #   withProgress(message = 'Initializing analysis', value = 0, {
    #     req(column())
    #     req(input$panel_size)
    #     req(sce())
    #     
    #     incProgress(detail = "Recomputing markers")
    #     
    #     ## Recompute the set of allowed genes (antibody applications might have changed) and get the markers
    #     allowed_genes(
    #       get_allowed_genes(input$select_aa, applications_parsed, sce())
    #     )
    #     
    #     ## Get the markers first time          
    #     fms(
    #       compute_fm(sce(), 
    #                  column(), 
    #                  pref_assay(),
    #                  allowed_genes()
    #       )
    #     )
    #     # 
    #     
    #     markers <- get_markers(fms(), 
    #                            input$panel_size, 
    #                            input$marker_strategy, 
    #                            sce())
    #     
    #     markers <- list(recommended_markers = input$bl_recommended,
    #          scratch_markers = input$bl_scratch,
    #          top_markers = input$bl_top)
    #     
    #     # SMH
    #     current_markers(set_current_markers_safely(markers, fms()))
    #     
    #     num_selected(length(current_markers()$top_markers))
    #   })
    #   
    #   update_analysis()
    # })
    
    observeEvent(input$add_cell_type_markers, {
      if(!is.null(input$cell_type_markers) && !is.null(fms())) {
        tmp <- get_cell_type_add_markers_reactable(fms()[[1]][[input$cell_type_markers]],
                                            unique(unlist(current_markers())))
        
        current_cell_type_marker_fm <<- tmp$fm
        marker_suggestions(tmp$fm)
        
        output$cell_type_marker_reactable <- renderReactable({ tmp$reactable })
        
        modal_add_marker <- modalDialog(
          reactableOutput('cell_type_marker_reactable'),
          title = paste("Select markers for", input$cell_type_markers),
          size = "m",
          easyClose = TRUE,        
          footer = tagList(
            modalButton("Cancel"),
            actionButton("add_select_markers", "Add selected")
          )
        )
        showModal(modal_add_marker)
      }
    })
    
    observeEvent(input$add_select_markers, {
      genes_to_add <- current_cell_type_marker_fm[selected_cell_type_markers(),]$Gene
      
      ## The following should really be turned into a function given
      ## frequently it is called
      cm <- current_markers()
      markers <- list(recommended_markers = cm$recommended_markers,
                      scratch_markers = input$bl_scratch,
                      top_markers = unique(c(genes_to_add, cm$top_markers)))
      current_markers(
        set_current_markers_safely(markers, fms())
      )
      
      num_markers_in_selected(length(current_markers()$top_markers))
      num_markers_in_scratch(length(current_markers()$scratch_markers))
      unique_cell_types <- unique(unlist(current_markers()$associated_cell_types))
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                unique_cell_types)
      
      removeModal()
    })
    
    ### MARKER SELECTION ###
    observeEvent(input$enter_marker, { # Manually add markers one by one
      
      if(!is.null(input$add_markers) && stringr::str_length(input$add_markers) > 1 && 
         (input$add_markers %in% allowed_genes()) && (! input$add_markers %in% 
            current_markers()$top_markers) && (! input$add_markers %in% 
                                               current_markers()$scratch_markers)) {
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
        
        num_markers_in_selected(length(current_markers()$top_markers))
        num_markers_in_scratch(length(current_markers()$scratch_markers))
        unique_cell_types <- unique(unlist(current_markers()$associated_cell_types))
        update_BL(current_markers(), num_markers_in_selected(),
                  num_markers_in_scratch(),
                  unique_cell_types)
        
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
      
      num_markers_in_selected(length(current_markers()$top_markers))
      num_markers_in_scratch(length(current_markers()$scratch_markers))
      unique_cell_types <- unique(unlist(current_markers()$associated_cell_types))
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                unique_cell_types)
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
      
      num_markers_in_selected(length(current_markers()$top_markers))
      num_markers_in_scratch(length(current_markers()$scratch_markers))
      unique_cell_types <- unique(unlist(current_markers()$associated_cell_types))
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                unique_cell_types)
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
              create_heatmap(sce()[,sce()$keep_for_analysis == "Yes"], current_markers(), col, display(), "Expression", pref_assay())
            )
          }
        } else { # Display heatmap for gene expression in a specific column
          
          for(col in columns) {
            if(col == display()) {
              heatmap(
                create_heatmap(sce()[,sce()$keep_for_analysis == "Yes"], current_markers(), display(), display(), input$heatmap_expression_norm, pref_assay())
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
      
      expression <- as.matrix(assay(sce()[,sce()$keep_for_analysis == "Yes"], pref_assay())[current_markers()$top_markers,])
      cmat <- cor(t(expression))
      
      # suggestions <- suggest_genes_to_remove(cmat, input$n_genes)
      
      suggestions(suggest_genes_to_remove(cmat, input$n_genes))
      
      showModal(suggestion_modal(suggestions = suggestions()))
    })
    
    observeEvent(input$remove_suggested, { # Remove suggested markers
      cm <- current_markers()
      
      if (!is.null(input$markers_to_remove)) {
        remove_marker <- c(input$markers_to_remove)
        num_markers_in_selected(length(current_markers()$top_markers))
        num_markers_in_scratch(length(current_markers()$scratch_markers))

        markers <- list(recommended_markers = cm$recommended_markers,
                   scratch_markers = unique(c(remove_marker, input$bl_scratch)),
                   top_markers = setdiff(cm$top_markers, remove_marker))       
        
        current_markers(
          set_current_markers_safely(markers, fms())
        )
        
        num_markers_in_selected(length(current_markers()$top_markers))
        num_markers_in_scratch(length(current_markers()$scratch_markers))
        unique_cell_types <- unique(unlist(current_markers()$associated_cell_types))
        update_BL(current_markers(), num_markers_in_selected(),
                  num_markers_in_scratch(),
                  unique_cell_types)
        
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
            compute_alternatives(input$input_gene, 
                                 sce()[,sce()$keep_for_analysis == "Yes"], 
                                 pref_assay(), input$number_correlations) |>
              drop_na()
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
      replacements <- as.vector(replacements()[n,]$Gene)
      alternative_marks(replacements)
      cm <- current_markers()
      markers <- list(recommended_markers = cm$recommended_markers,
                      scratch_markers = input$bl_scratch,
                      top_markers = unique(c(replacements, setdiff(cm$top_markers, input$bl_scratch))))
      
      # SMH
      current_markers(
        set_current_markers_safely(markers, fms())
      )
      
      num_markers_in_selected(length(current_markers()$top_markers))
      num_markers_in_scratch(length(current_markers()$scratch_markers))
      
      unique_cell_types <- unique(unlist(current_markers()$associated_cell_types))
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                unique_cell_types, adding_alternative = F)
      
      showNotification("Marker(s) added successfully.",
                       duration = 3)
      
    })
    
    ### UPDATE SELECTED MARKERS ###
    update_BL <- function(markers, top_size, scratch_size, unique_cell_types,
                          adding_alternative = F) {
      

      # unique_cell_types <- sort(unique(markers$associated_cell_types))
      # n_cell_types <- length(unique_cell_types)
      # palette <<- sample(cell_type_colors)[seq_len(n_cell_types)]
      set.seed(12345L)
      
      unique_cell_types <- sort(unique_cell_types)
      
      palette_to_use <- full_palette[1:length(unique_cell_types)]
      names(palette_to_use) <- unique_cell_types
      
      cytosel_palette(palette_to_use)
      
      # markers$top_markers <- sapply(markers$top_markers, function(m) paste(icon("calendar"), m))
      
      labels_top <- lapply(markers$top_markers, 
                           function(x) div(x, map_gene_name_to_antibody_icon(x, markers), style=paste('padding: 3px; color:', 
                                                                                                      set_text_colour_based_on_background(cytosel_palette()[ markers$associated_cell_types[x]]), '; background-color:', 
                                                                                                      cytosel_palette()[ markers$associated_cell_types[x] ])))
      labels_scratch <- lapply(markers$scratch_markers, 
                               function(x) div(x, map_gene_name_to_antibody_icon(x, markers), style=paste('padding: 3px; color:', 
                                                                                                          set_text_colour_based_on_background(cytosel_palette()[ markers$associated_cell_types[x]]), '; background-color:', 
     
                                                                                   cytosel_palette()[ markers$associated_cell_types[x] ])))
 
      output$legend <- renderPlot(cowplot::ggdraw(get_legend(cytosel_palette())))
     
      
      output$BL <- renderUI({
        bucket_list(
          header = "Marker selection",
          orientation = "horizontal",
          group_name = "bucket_list_group",
          add_rank_list(
            text = paste0("Scatch space (", scratch_size, ")"),
            labels = labels_scratch,
            input_id = "bl_scratch",
            options = c(multiDrag = TRUE),
            class = c("default-sortable", "cytocellbl")
          ),
          add_rank_list(
            text = paste0("Selected markers (", top_size, ")"),
            labels = labels_top,
            input_id = "bl_top",
            options = c(multiDrag = TRUE),
            class = c("default-sortable", "cytocellbl")
          )
        )
      })
      
    }
    
    
    ### UPDATE ANALYSIS ###
    update_analysis <- function() {
      
      withProgress(message = 'Processing data', value = 0, {

        
        ## Re-set the set of allowed genes (these may have changed if a different
        ## antibody application is selected)
        allowed_genes(
          get_allowed_genes(input$select_aa, applications_parsed, 
                            sce()[,sce()$keep_for_analysis == "Yes"])
        )
        
        ## Set that these genes can be selected
        update_autocomplete_input(session, "add_markers",
                                  options = allowed_genes())
        
        update_autocomplete_input(session, "input_gene",
                                  options = allowed_genes())
        
        
        num_markers_in_selected(length(current_markers()$top_markers))
        num_markers_in_scratch(length(current_markers()$scratch_markers))
        
        
        unique_cell_types <- unique(unlist(current_markers()$associated_cell_types))
        
        update_BL(current_markers(), num_markers_in_selected(),
                  num_markers_in_scratch(),
                  unique_cell_types)
        
        # Update UMAP
        incProgress(detail = "Computing UMAP")
        umap_all(get_umap(sce()[,sce()$keep_for_analysis == "Yes"],
                          column(), pref_assay()))
        umap_top(get_umap(sce()[,sce()$keep_for_analysis == "Yes"][current_markers()$top_markers,], column(), pref_assay()))
        
        
        # Update heatmap
        incProgress(detail = "Drawing heatmap")
        
        columns <- column()
        toggle(id = "norm", condition = display() != "Marker-marker correlation")
        
        if(display() == "Marker-marker correlation") {
          for(col in columns) {
            heatmap(
              create_heatmap(sce()[,sce()$keep_for_analysis == "Yes"], 
                             current_markers(), col, 
                             display(), "Expression", pref_assay()))
            
          }
        } else {
          for(col in columns) {
            if(col == display()) {
              heatmap(
                create_heatmap(sce()[,sce()$keep_for_analysis == "Yes"], 
                               current_markers(), display(), 
                               display(), input$heatmap_expression_norm, 
                               pref_assay())
              )
              
              
              break
            }
          }
        }
        cells_per_type(table(colData(sce()[,
                            sce()$keep_for_analysis == "Yes"])[[column()]]))
        
        # Update metrics
        incProgress(detail = "Computing panel score")
        scores <- get_scores(sce()[,sce()$keep_for_analysis == "Yes"], 
                             column(), current_markers()$top_markers, pref_assay())
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
  

  
  shinyApp(ui, server, ...)
}


