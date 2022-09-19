# utils::globalVariables(c("Cell type", "Symbol", "r", "score", "what"), "cytosel")

# palette <- NULL

utils::globalVariables(c("palette_to_use", "full_palette"), "cytosel")

ggplot2::theme_set(cowplot::theme_cowplot())

options(shiny.maxRequestSize = 1000 * 200 * 1024 ^ 2, warn=-1)

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
#' @importFrom dplyr desc mutate_if
#' @importFrom bsplus use_bs_popover shinyInput_label_embed shiny_iconlink bs_embed_tooltip use_bs_tooltip
#' @importFrom shinyjs useShinyjs hidden toggle
#' @importFrom grDevices dev.off pdf
#' @importFrom zip zip
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom SingleCellExperiment reducedDimNames reducedDims
#' @importFrom plotly plot_ly plotlyOutput renderPlotly layout
#' @importFrom parallelly availableCores
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stringr str_split_fixed
#' @importFrom shinydashboard box dashboardBody dashboardHeader dashboardSidebar
#' dashboardPage menuItem sidebarMenu sidebarMenuOutput tabItem tabItems
#' valueBoxOutput renderMenu updateTabItems
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom yaml read_yaml
#' @export
#' 
#' @param ... Additional arguments
cytosel <- function(...) {
  
  antibody_info <- dplyr::rename(cytosel_data$antibody_info, Symbol = `Gene Name (Upper)`)
  antibody_info <- tidyr::drop_na(antibody_info)
  
  options(MulticoreParam=quote(MulticoreParam(workers=availableCores())))

  
  applications_parsed <- get_antibody_applications(antibody_info, 
                                                   'Symbol', 'Listed Applications')
  
  grch38 <- cytosel_data$grch38
  
  full_palette <- create_global_colour_palette()
  
  ui <- tagList(
    includeCSS(system.file(file.path("www", "cytosel.css"),
                          package = "cytosel")),
    tags$head(
      tags$script(HTML("window.onbeforeunload = function() {return 'Your changes will be lost!';};"))
    ),
    tags$head(tags$style(".modal-dialog{ width:750px}")),
    tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),

    # Use packages
    useShinyalert(force = TRUE),
    use_bs_tooltip(),
    useShinyjs(),
    
    dashboardPage(
    
    # Title
    dashboardHeader(title = "cytosel"
                    ),
    
    dashboardSidebar(
      sidebarMenu(id = "tabs",
        # width = 300,
        menuItem("Get Started", tabName = "inputs", icon = icon("gear")),
        sidebarMenuOutput(outputId = 'output_menu'),
        br(),
        # div(style="display:inline-block;width:110%;text-align: center;",
        #     actionButton("start_analysis", "Go!", icon=icon("play", style = "color:black;"))),
        column(10, offset = 0, align = "center",
                        actionButton("start_analysis", "Run analysis!", icon=icon("play", style = "color:black;"))),
        column(10, offset = 0, align = "center",
            downloadButton("downloadData", "Save panel", style = "color:black;
                           margin-top: 15px;"))
      )
    ),
    dashboardBody(
      
      # https://stackoverflow.com/questions/52198452/how-to-change-the-background-color-of-the-shiny-dashboard-body
      tags$head(tags$style(HTML('
                                /* body */
                                .content-wrapper, .right-side {
                                background-color: #FFFFFF;
                                }
                                
                                '))),
      
      
    # bsPopover(id = 'q1', 'Upload an Input scRNA-seq file', content = 'Input scRNA-seq file',
    #             trigger = 'click', options = list(container = "body")),
      
    # Tabs
    tabItems(
      tabItem("inputs",
              fluidRow(column(5,
        fileInput("input_scrnaseq",
                  label = p(
                    'Input scRNA-seq',
                  ), accept = c(".rds")) %>%
          shinyInput_label_embed(
            icon("circle-info") %>%
              bs_embed_tooltip(title = get_tooltip('input_scrnaseq'),
                               placement = "right", html = "true")
          ))),
        selectInput("coldata_column", "Cell category to evaluate:", NULL, multiple=FALSE) %>%
          shinyInput_label_embed(
            icon("circle-info") %>%
              bs_embed_tooltip(title = get_tooltip('coldata_column'),
                      placement = "right", html = "true")
          ),
        # add padding space between elements
        fluidRow(column(4, bsCollapse(id = "advanced_collapse",
                       bsCollapsePanel(title = HTML(paste0(
                         "Advanced settings", tags$span(icon("sort-down",
                                                       style = "position:right; margin-left: 4px; margin-top: -4px;")))), style = "info",
                                   textOutput("selected_assay"),
                                   br(),
                                   actionButton("show_cat_table", "Category subsetting",
                                    align = "center",
                                    width = NULL),
                   br(),
                   br(),
                   numericInput("panel_size", "Targeted panel size:", 32, 
                                min = 1, max = 200, step = NA, width = NULL) %>%
                                shinyInput_label_embed(
                                icon("circle-info") %>%
                                bs_embed_tooltip(title = get_tooltip('panel_size'),
                                placement = "right")),
                   radioButtons("marker_strategy", label = "Marker selection strategy",
                                choices = list("Cell type based"="fm", "Cell type free (geneBasis)" = "geneBasis"),
                                selected="fm"),
                   checkboxInput("subsample_sce", "Subsample cells", value = TRUE) %>%
                     shinyInput_label_embed(
                     icon("circle-info") %>%
                     bs_embed_tooltip(title = get_tooltip('subsample_sce'),
                     placement = "right")),
                   hidden(div(id = "precomputed",
                   checkboxInput("precomputed_dim", "Use precomputed UMAP", value = F))),
                   selectInput("select_aa", "Antibody applications:", 
                              applications_parsed$unique_applications, multiple=TRUE) %>%
                               shinyInput_label_embed(
                               icon("circle-info") %>%
                               bs_embed_tooltip(title = "Placeholder",
                               placement = "right", html = "true"))
                   )))),
        fluidRow(column(5,
                        fileInput("read_back_analysis",
                                  label = p(
                                    'Upload previous analysis',
                                  ), accept = c(".yml")) %>%
                          shinyInput_label_embed(
                            icon("circle-info") %>%
                              bs_embed_tooltip(title = get_tooltip('reupload'),
                                               placement = "right", html = "true")
                          )))),
      tabItem("marker_selection",
                 # icon = icon("list"),
                 br(),
                 fluidRow(column(12, actionButton("markers_change_modal", "Add markers to panel"))),
                 hr(),
                 hidden(div(id = "marker_display",
                   tags$span(icon("circle-info")
                           %>%
                   bs_embed_tooltip(title = get_tooltip('marker_display'),
                                    placement = "right", html = TRUE)))),
                 hidden(div(id = "marker_visualization",
                            tags$span(icon("circle-info")
                                      %>%
                                        bs_embed_tooltip(title = get_tooltip('marker_visualization'),
                                                         placement = "right", html = TRUE)))),
                fluidRow(column(4, align = "center", div(style="display: inline-block; font-size: 15px", 
                                  htmlOutput("scratch_marker_counts"))),
                         column(4, align = "center", div(style="display: inline-block; font-size: 15px", 
                                  htmlOutput("selected_marker_counts")))),
                 fluidRow(column(8, uiOutput("BL"),
                                 style = "margin-bottom:0px;"
                                 ),
                          column(3, plotOutput("legend", width = "250px"),
                                 style = "margin-top:-25px;"))),
        
      tabItem("UMAP",
                 # icon = icon("globe"),
                 br(),
                 helpText(get_tooltip('umap')),
              br(),
              br(),
                 fluidRow(column(6, plotlyOutput("all_plot", width="500px", height="350px")),
                          column(6, plotlyOutput("top_plot", width="500px", height="350px")))
          ),
        
      tabItem("Heatmap",
                 # icon = icon("table"),
                 br(),
                 splitLayout(cellWidths = c(320, 280),
                             div(selectInput("display_options", 
                              "Display expression or gene correlation", 
                              choices = c("Marker-marker correlation"), width = "86%") %>%
                                 shinyInput_label_embed(icon("circle-info") %>%
                                bs_embed_tooltip(title = get_tooltip('heatmap_display_options'),
                                  placement = "right", html = "true"))
                              ),
                             hidden(div(id = "norm",
                                selectInput("heatmap_expression_norm", 
                                            "Heatmap expression normalization", 
                                            choices = c("Expression", "z-score"), width = "89%") %>%
                                  shinyInput_label_embed(icon("circle-info") %>%
                                  bs_embed_tooltip(title = get_tooltip('heatmap_expression_norm'),
                                      placement = "right"))
                                ))),
              splitLayout(cellWidths = c(190, 200),
                          div(numericInput("n_genes", "Number of genes to remove", 
                                           value = 10, min = 1, width = "50%")),
                          div(style = "margin-top:25px;", actionButton("suggest_gene_removal", "View suggestions"))),
                 plotlyOutput("heatmap", height="600px"),
                 hr(),
                 tags$span(icon("circle-info") %>%
                             bs_embed_tooltip(title = get_tooltip('gene_removal'),
                                              placement = "right"))
          ),

      tabItem("Metrics",
                 # icon = icon("ruler"),
                 br(),
                 helpText(get_tooltip('metrics')),
                 tags$span(icon("circle-info") %>%
                             bs_embed_tooltip(title = get_tooltip('metrics_explanation'),
                                              placement = "right", html = TRUE)),
                 textOutput("cells_per_category"),
                 plotlyOutput("metric_plot")
          ),

      tabItem("alternative_markers",
                 # icon = icon("exchange-alt"),
                 br(),
                 helpText(get_tooltip("alternative_markers")),
                 splitLayout(cellWidths = c(180, 240),
                             div(autocomplete_input("input_gene", "Input gene", options=c(), width = "100%")),
                             div(numericInput("number_correlations", "Number of alternative markers", 
                                              value = 10, min = 1, width = "35%")),
                             hidden(div(id = "send", actionButton("send_markers", 
                                    "Send markers to selection panel")))),
                 actionButton("enter_gene", "Enter"),
                 br(),
                 br(),
                 DTOutput("alternative_markers")
                 # hidden(div(id = "send", actionButton("send_markers", "Send markers to selection panel"))),
                 # br()
          ),

      tabItem("antibody_explorer",
                 # icon = icon("wpexplorer"),
                 br(),
                 reactableOutput("antibody_table")
          )
      )
    )
  )
  )
  
  server <- function(input, output, session) {
    
    ### REACTIVE VARIABLES ###
    plots <- reactiveValues() # Save plots for download
    
    umap_top <- reactiveVal()
    umap_all <- reactiveVal()
    heatmap <- reactiveVal()
    metrics <- reactiveVal()
    previous_metrics <- reactiveVal()
    current_metrics <- reactiveVal()
    
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
    
    cell_min_threshold <- reactiveVal()
    
    use_precomputed_umap <- reactiveVal(FALSE)
    possible_umap_dims <- reactiveVal()
    umap_precomputed_col <- reactiveVal(NULL)
    first_render_outputs <- reactiveVal(FALSE)
    df_antibody <- reactiveVal()
    
    view_advanced_settings <- reactiveVal(FALSE)
    
    reupload_analysis <- reactiveVal(FALSE)
    markers_reupload <- reactiveVal()
    yaml <- reactiveVal()
    
    # addPopover(session, "q1", "Upload an Input scRNA-seq file", content = 'Input scRNA-seq file',
    #            trigger = 'click')
    
    ### UPLOAD FILE ###
    observeEvent(input$input_scrnaseq, {
      #if(isTruthy(methods::is(obj, 'SingleCellExperiment')) || isTruthy(methods::is(obj, 'Seurat'))) {
      
      updateCheckboxInput(inputId = "precomputed_dim", value = F)
      use_precomputed_umap(FALSE)
      umap_precomputed_col(NULL)
      input_sce <- read_input_scrnaseq(input$input_scrnaseq$datapath)
      input_sce <- detect_assay_and_create_logcounts(input_sce)
      input_sce <- parse_gene_names(input_sce, grch38)
      sce(input_sce)
      input_assays <- c(names(assays(sce())))
      
      if (!isTruthy(input$min_category_count)) {
        updateNumericInput(session, "min_category_count", 2)
        cell_min_threshold(2)
      } else {
        cell_min_threshold(input$min_category_count)
      }
    
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
      
      if (!isTruthy(input$coldata_column)) {
        column(colnames(colData(sce()))[1])
      } else {
        column(input$coldata_column)
      }
      
      possible_umap_dims(detect_umap_dims_in_sce(sce()))
      
      toggle(id = "precomputed", condition = length(possible_umap_dims()) > 0)

    })
    
    observeEvent(input$coldata_column, {
      req(input$input_scrnaseq)
      column(input$coldata_column)
      
      if (!isTruthy(reupload_analysis())) {
        specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      }
      })
    
    observeEvent(input$add_selected_to_analysis, {
      req(input$input_scrnaseq)
      req(input$coldata_column)
      specific_cell_types_selected(input$user_selected_cells)
      
      removeModal()
      
      if (length(specific_cell_types_selected()) > 0) {
        showNotification("Setting subset to select cell types",
                         duration = 2)
      } else {
        showNotification("Empty subset selection, defaulting to using all cell types in analysis.",
                         duration = 2)
        specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      }
      
    })
    
    observeEvent(input$user_selected_cells, {
      req(input$input_scrnaseq)
      req(input$coldata_column)
      specific_cell_types_selected(input$user_selected_cells)
    })
    
    observeEvent(input$min_category_count, {
      cell_min_threshold(input$min_category_count)
    })
    
    observeEvent(input$show_cat_table, {
    req(sce())
    req(pref_assay())
    req(input$coldata_column)
    # req(cell_min_threshold())
  
    cell_category_table <- create_table_of_hetero_cat(sce(),
                                                      input$coldata_column)

      if (!isTruthy(cell_min_threshold())) {
        cell_min_threshold(2)
      }
  
      if (!isTruthy(specific_cell_types_selected())) {
        specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      }
  
    
    if (nrow(cell_category_table) > 100) {
      
      too_large_to_show_modal(input$coldata_column)
    } else {
      summary_cat_tally(cell_category_table)
      
      output$cell_cat_table <- DT::renderDataTable(
        summary_cat_tally()
      )
      
      showModal(cell_cat_modal(cell_min_threshold(), unique(sce()[[input$coldata_column]]),
                               specific_cell_types_selected()))
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
    
    ### bring up modal to add markers ###
    observeEvent(input$markers_change_modal, {
      req(sce())
      req(current_markers())
      req(allowed_genes())
      req(fms())
      
      # updateSelectInput(
      #   session = session,
      #   inputId = "cell_type_markers",
      #   choices = names(fms()[[1]])
      # )
      # 
      # updateSelectizeInput(
      #   session = session,
      #   inputId = "add_markers",
      #   choices = c("", allowed_genes())
      # )
    
      showModal(markers_add_modal(allowed_genes(), names(fms()[[1]])))
    })
    
    ### ANTIBODY EXPLORER ###
    output$antibody_table <- renderReactable({
      req(current_markers())
      req(df_antibody())
  
      reactable(df_antibody(),
                searchable = TRUE,
                filterable = TRUE,
                sortable = TRUE)
    })
    
    ### PLOTS ###
    output$all_plot <- renderPlotly({
      req(umap_all)
      req(column)
      req(plots$all_plot)
      
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

      plots$all_plot
    })
    # },
    # width=350,
    # height=300)
    
    output$top_plot <- renderPlotly({
      req(umap_top)
      req(column)
      req(plots$top_plot)
      
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
      
      plots$top_plot
      
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
      req(plots$metric_plot)
     
      plots$metric_plot
      
    })
    
    ### ANALYSIS ###
    observeEvent(input$start_analysis, {
      req(input$panel_size)
      req(input$coldata_column)
      req(sce())
      
      if (isTruthy(current_metrics())) {
        previous_metrics(current_metrics() |> mutate(Run = "Previous"))
      }
      
      if (!isTruthy(specific_cell_types_selected())) {
          specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      }
        
        if (!isTruthy(cell_min_threshold())) {
          cell_min_threshold(2)
        }
      
      if (cell_min_threshold() < 2) {
        proceed_with_analysis(FALSE)
        threshold_too_low_modal()
      } else {
        proceed_with_analysis(TRUE)
      }
      
      # if the user never opened up the tally table, automatically set all cell types
      # for analysis
      
      withProgress(message = 'Initializing analysis', value = 0, {
        incProgress(detail = "Checking data")
        
        cell_cat_summary <- create_table_of_hetero_cat(sce(),
                                                       input$coldata_column)
        
        cell_types_high_enough(remove_cell_types_by_min_counts(cell_cat_summary,
                                                               sce(), 
                                                               input$coldata_column, 
                                                               cell_min_threshold()))
        
        cell_types_to_keep(intersect(specific_cell_types_selected(),
                                     cell_types_high_enough()))
        
        sce(remove_null_and_va_from_cell_cat(sce(), input$coldata_column))
        
        showNotification("Ignoring any cells with input category set to NA or null",
                         type = 'message',
                         duration = 4)
        
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
          cell_type_ignored_modal(cell_min_threshold(),
                                  cell_types_excluded())
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
          
          high_cell_number_warning(ncol(sce()[,sce()$keep_for_analysis == "Yes"]), 10000)
          
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
            
            if(!is.null(input$bl_top)) {
              ## We get here if input$bl_top exists, ie if this
              ## is an analysis refresh
              ## in this case we set the markers to their existing values
              markers <- list(recommended_markers = input$bl_recommended,
                              scratch_markers = input$bl_scratch,
                              top_markers = input$bl_top)
            } else {
              ## We compute the set of markers for the first time
              
              if (isTruthy(markers_reupload())) {
                markers <- list(recommended_markers = markers_reupload()$top_markers,
                                top_markers = markers_reupload()$top_markers,
                                scratch_markers = markers_reupload()$scratch_markers)
              } else {
                markers <- get_markers(fms(), 
                                       # Adding 10 to make sure panel size is approximate 
                                       # since a) the same marker is selected multiple times and
                                       # b) excess markers are removed
                                       input$panel_size, 
                                       input$marker_strategy, 
                                       sce()[,sce()$keep_for_analysis == "Yes"],
                                       allowed_genes())
                
                if(length(markers$recommended_markers) < input$panel_size){
                  showNotification("Cytosel found genes that are good markers for multiple cell types. This will result in a smaller panel size than requested.
                                 You may manually add additional markers.",
                                   type = 'message',
                                   duration = NULL)
                }
                ## Forgotten what this is for
                markers$scratch_markers <- scratch_markers_to_keep
              }
            }
            
            # SMH
            current_markers(set_current_markers_safely(markers, fms()))
          
            num_markers_in_selected(length(current_markers()$top_markers))
            num_markers_in_scratch(length(current_markers()$scratch_markers))
            cells_per_type(table(colData(
              sce()[,sce()$keep_for_analysis == "Yes"])[[column()]]))
            
            df_antibody(dplyr::filter(antibody_info, Symbol %in% 
                                        current_markers()$top_markers))
            
            update_analysis()
            
            
          } else {
            unique_element_modal(col)
          }
        
          if (!isTruthy(first_render_outputs()) & isTruthy(current_markers())) {
            
            output$output_menu <- renderMenu(expr = {
              sidebarMenu(
                menuItem("Marker selection", tabName = "marker_selection", 
                         icon = icon("barcode")),
                menuItem("UMAP", tabName = "UMAP", 
                         icon = icon("arrows-alt")),
                menuItem("Heatmap", tabName = "Heatmap", 
                         icon = icon("th-large")),
                menuItem("Metrics", tabName = "Metrics", 
                         icon = icon("line-chart")),
                menuItem("Alternative Markers", tabName = "alternative_markers", 
                         icon = icon("arrows-h")),
                menuItem("Antibody Explorer", tabName = "antibody_explorer", 
                         icon = icon("list-alt"))
                )
            })
            
            first_render_outputs(TRUE)
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
    observeEvent(input$bl_top, {
      req(sce())
      req(current_markers())
      
      current_markers(list(recommended_markers = current_markers()$recommended_markers,
                      scratch_markers = input$bl_scratch,
                      top_markers = input$bl_top,
                      associated_cell_types = current_markers()$associated_cell_types))
      
      # current_markers(set_current_markers_safely(markers, fms()))
      num_markers_in_selected(length(current_markers()$top_markers))
      
      output$selected_marker_counts<- renderText({paste("<B>",
                                      "Selected Markers:", num_markers_in_selected(), "</B>")})
      
    })
    
    observeEvent(input$bl_scratch, {
      req(sce())
      req(current_markers())
      
      current_markers(list(recommended_markers = current_markers()$recommended_markers,
                           scratch_markers = input$bl_scratch,
                           top_markers = input$bl_top,
                           associated_cell_types = current_markers()$associated_cell_types))
      
      # current_markers(set_current_markers_safely(markers, fms()))
      num_markers_in_scratch(length(current_markers()$scratch_markers))
      
      output$scratch_marker_counts <- renderText({paste("<B>", 
                                                        "Scratch Markers:", num_markers_in_scratch(), "</B>")})
      
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
      
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                names(fms()[[1]]))
      
      removeModal()
    })
    
    ### MARKER SELECTION ###
    observeEvent(input$enter_marker, { # Manually add markers one by one
      req(input$add_markers)
      
      if(!is.null(input$add_markers) &&
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
        
        update_BL(current_markers(), num_markers_in_selected(),
                  num_markers_in_scratch(),
                  names(fms()[[1]]))
        
        removeModal()
        
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
      
      
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                names(fms()[[1]]))
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
      
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                names(fms()[[1]]))
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
      req(current_markers())
      
      expression <- as.matrix(assay(sce()[,sce()$keep_for_analysis == "Yes"], pref_assay())[current_markers()$top_markers,])
      cmat <- cor(t(expression))
      
      # suggestions <- suggest_genes_to_remove(cmat, input$n_genes)
      
      suggestions(suggest_genes_to_remove(cmat, input$n_genes))
      
      showModal(suggestion_modal(suggestions = suggestions(),
                                 possible_removal = current_markers()$top_markers))
    })
    
    observeEvent(input$remove_suggested, { # Remove suggested markers
      req(current_markers())
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
        
        update_BL(current_markers(), num_markers_in_selected(),
                               num_markers_in_scratch(),
                               names(fms()[[1]]))
        
        removeModal()
        
        # updateTabItems(session, "tabs", "inputs")
        updateTabsetPanel(session, "tabs", "marker_selection")
        
      } else {
        showModal(suggestion_modal(failed = TRUE, suggestions()))
      }
    
    })
    
    
    ### ALTERNATIVE MARKERS ###
    observeEvent(input$enter_gene, { # Compute alternative markers
      req(input$number_correlations)
      req(sce())
      
      if(!is.null(input$input_gene) && stringr::str_length(input$input_gene) > 1 && (input$input_gene %in% 
                                    rownames(sce()[,sce()$keep_for_analysis == "Yes"])) &&
         (input$input_gene %in% allowed_genes())) {
        
        withProgress(message = 'Updating analysis', value = 0, {
          incProgress(6, detail = "Computing alternatives")
          
          # Make this sampling dependent on the input sample argument
          replacements(
            compute_alternatives(input$input_gene, 
                                 sce()[,sce()$keep_for_analysis == "Yes"], 
                                 pref_assay(), input$number_correlations,
                                 allowed_genes()) |>
              drop_na()
          )
          
          output$alternative_markers <- renderDT(replacements(), server = TRUE)
        })
        
      } else if(!(input$input_gene %in% rownames(sce()[,sce()$keep_for_analysis == "Yes"])) |
                !(input$input_gene %in% allowed_genes())) {
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
  
      
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                names(fms()[[1]]))
      
      showNotification("Marker(s) added successfully.",
                       duration = 3)
      
      updateTabsetPanel(session, "tabs", "marker_selection")
      
      
    })
    
    observeEvent(input$precomputed_dim, {
      req(sce())
      req(possible_umap_dims())
      
      if (isTruthy(input$precomputed_dim)) {
        showModal(pre_computed_umap_modal(possible_umap_dims()))
      }
    })
    
    observeEvent(input$select_precomputed_umap, {
      req(sce())
      req(input$possible_precomputed_dims)
      
      if (isTruthy(input$precomputed_dim)) {
        use_precomputed_umap(TRUE)
        umap_precomputed_col(input$possible_precomputed_dims)
        removeModal()
      }
      
    })
    
    ## Re-upload previous analysis
    
    observeEvent(input$read_back_analysis, {
      req(sce())
      
      yaml_back <- read_back_in_saved_yaml(input$read_back_analysis$datapath)
      
      if (isTruthy(yaml_back$`Pre-computed UMAP`) &
          length(reducedDimNames(sce())[grepl("UMAP|umap|Umap|uMap|uMAP",
          reducedDimNames(sce()))]) < 1) {
        reupload_failed_modal()
        reupload_analysis(FALSE)
      } else if (yaml_back$`Number of columns (cells)` == ncol(sce()) &
          yaml_back$`Number of rows (features)` == nrow(sce()) &
          yaml_back$`Heterogeneity source` %in% colnames(colData(sce()))) {
        
        updateSelectInput(session, "coldata_column", choices = colnames(colData(sce())),
                          selected = yaml_back$`Heterogeneity source`)
        pref_assay(yaml_back$`Assay used`)
        updateNumericInput(session, "panel_size", value = yaml_back$`Target panel size`)
        cell_min_threshold(yaml_back$`Min Cell Category cutoff`)
        updateNumericInput(session, "min_category_count", value = yaml_back$`Min Cell Category cutoff`)
        updateCheckboxInput(session, "subsample_sce", value = yaml_back$`Subsampling Used`)
        updateRadioButtons(session, "marker_strategy", selected = yaml_back$`Marker strategy`)
        updateSelectInput(session, "select_aa", selected = yaml_back$`Antibody applications`)
        specific_cell_types_selected(yaml_back$`User selected cells`)
        markers_reupload(list(top_markers = yaml_back$`Selected marker panel`,
                              scratch_markers = yaml_back$`Scratch marker panel`))
        yaml_length <- sum(c(length(markers_reupload()$top_markers),
                             length(markers_reupload()$scratch_markers)))
        if (yaml_length != yaml_back$`Target panel size`) {
          updateNumericInput(session, "panel_size", value = yaml_length)
        }
        updateCheckboxInput(session, inputId = "precomputed_dim",
                            value = yaml_back$`Pre-computed UMAP`)
        use_precomputed_umap(yaml_back$`Pre-computed UMAP`)
        
        if (isTruthy(use_precomputed_umap())) {
          possible_umap_dims(detect_umap_dims_in_sce(sce()))
          showModal(pre_computed_umap_modal(possible_umap_dims()))
        }
        
        reupload_analysis(TRUE)
      } else {
        reupload_failed_modal()
        reupload_analysis(FALSE)
      }
      
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
          header = NULL,
          orientation = "horizontal",
          group_name = "bucket_list_group",
          add_rank_list(
            text = NULL,
            labels = labels_scratch,
            input_id = "bl_scratch",
            options = c(multiDrag = TRUE),
            class = c("default-sortable", "cytocellbl")
          ),
          add_rank_list(
            labels = labels_top,
            text = NULL,
            input_id = "bl_top",
            options = c(multiDrag = TRUE),
            class = c("default-sortable", "cytocellbl")
          )
        )
      })
      
    }
    
    
    ### UPDATE ANALYSIS ###
    update_analysis <- function() {
      
      withProgress(message = 'Updating analysis', value = 0, {

        
        ## Re-set the set of allowed genes (these may have changed if a different
        ## antibody application is selected)
        allowed_genes(
          get_allowed_genes(input$select_aa, applications_parsed, 
                            sce()[,sce()$keep_for_analysis == "Yes"])
        )
        
        ## Set that these genes can be selected
        update_autocomplete_input(session, "input_gene",
                                  options = allowed_genes())
        
        num_markers_in_selected(length(current_markers()$top_markers))
        num_markers_in_scratch(length(current_markers()$scratch_markers))
        
        update_BL(current_markers(), num_markers_in_selected(),
                  num_markers_in_scratch(),
                  names(fms()[[1]]))
        
        # Update UMAP
        incProgress(detail = "Computing & creating UMAP plots")
        umap_all(get_umap(sce()[,sce()$keep_for_analysis == "Yes"],
                          column(), pref_assay(), 
                          use_precomputed_umap(),
                          umap_precomputed_col(), F, num_markers_in_selected()))
        umap_top(get_umap(sce()[,sce()$keep_for_analysis == "Yes"][current_markers()$top_markers,], 
                          column(), pref_assay(), use_precomputed_umap(),
                          umap_precomputed_col(),
                          T, num_markers_in_selected()))
        
        plots$all_plot <- plot_ly(umap_all(), x=~UMAP1, y=~UMAP2, color=~get(columns[1]), text=~get(columns[1]), 
                                 type='scatter', hoverinfo="text", colors=cytosel_palette()) %>% 
          layout(title = "UMAP all genes", showlegend = F)
        
        plots$top_plot <- plot_ly(umap_top(), x=~UMAP1, y=~UMAP2, color=~get(columns[1]), text=~get(columns[1]), 
                                 type='scatter', hoverinfo="text", colors=cytosel_palette()) %>% 
          layout(title = "UMAP selected markers")
        
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
        
        print(metrics())
        
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
        
        current_metrics(m |> mutate(Run = "Current", cat = str_split_fixed(what, " ", 2)[,1]))
        
        if (isTruthy(previous_metrics())) {
          all_metrics <- rbind(previous_metrics(), current_metrics()) 
        } else {
          all_metrics <- current_metrics()
        }
        
        plots$metric_plot <- plot_ly(all_metrics, x = ~score, y = ~cat, 
                                     color = ~Run,
                                     type='box', hoverinfo = 'none') %>% 
          layout(boxmode = "group",
                 xaxis = list(title="Score"),
                 yaxis = list(title="Source"))
        
        
        # Show help text popover
        shinyjs::show(id = "marker_visualization")
        shinyjs::show(id = "marker_display")
        
        
      })
    }
    
    
    ### SAVE PANEL ###
    output$downloadData <- downloadHandler(
      filename <- paste0("Cytosel-Panel-", Sys.Date(), ".zip"),
      # reactive_vals <- c(current_markers(), heatmap(), pref_assay(),
      #                    column(), sce())
      # 
      # truthiness <- sapply(reactive_vals, FUN = function(x) isTruthy(x))
      
      content = function(fname) {
        download_data(fname, current_markers(), plots, heatmap(),
                      input_file = input$input_scrnaseq$datapath,
                      assay_used = pref_assay(),
                      het_source = column(),
                      panel_size = input$panel_size,
                      cell_cutoff_value = as.integer(cell_min_threshold()),
                      subsample = input$subsample_sce,
                      antibody_table = df_antibody(),
                      marker_strat = input$marker_strategy,
                      antibody_apps = input$select_aa,
                      selected_cell_types = input$user_selected_cells,
                      precomputed_umap_used = input$precomputed_dim,
                      num_cells = ncol(sce()),
                      num_genes = nrow(sce()))
      },
      contentType = "application/zip"
    )
    
  }
  

  
  shinyApp(ui, server, ...)
}


