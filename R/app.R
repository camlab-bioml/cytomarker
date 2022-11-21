
# utils::globalVariables(c("Cell type", "Symbol", "r", "score", "what"), "cytosel")

# palette <- NULL

### user selects pre-curated dataset ####

curated_dataset_names <- c("Seurat PBMC", "Baron et al. Pancreas")

preview_info <- c(paste("<b>", "<a href='https://satijalab.org/seurat/articles/pbmc3k_tutorial.html' target='_blank' >Seurat PBMC tutorial dataset</a>",
                        "</b>", "<br/>", "<b>", "Cells:", "</b>", "2638", 
                         "<br/>", "<b>", "Genes:", "</b>", "13,714", "<br/>", 
                        "<b>", "Cell category of interest:", "</b>", "ident",
                        "<br/>", "<b>", "Cell types in category:", "</b>", "<br/>", "Memory CD4 T, B, CD14+ Mono and 6 others"),
                  paste("<b>", "<a href='https://www.sciencedirect.com/science/article/pii/S2405471216302666' target='_blank' >Human pancreas, Baron et al. (2016)</a>",
                  "</b>",
                        "<br/>", "<b>", "Cells:", "</b>", "2069", 
                       "<br/>", "<b>", "Genes:", "</b>", "20,125", "<br/>", 
                       "<b>", "Cell category of interest:", "</b>", "label",
                       "<br/>", "<b>", "Cell types in category:", "</b>", "<br/>", "acinar, delta, beta and 5 others")) |>
  set_names(curated_dataset_names)

dataset_selections <- c("seurat_pbmc.rds", "baron_pancreas_ref.rds") |>
  set_names(curated_dataset_names)

default_celltype_curated <- c("ident", "label") |>
  set_names(curated_dataset_names)

# can use to remove any improperly downloaded data sets (shouldn't be necessary with refresh token)
# 
for (i in dataset_selections) {
  if (file.exists(file.path(tempdir(), "/", i))) {
    command <- paste('rm ', tempdir(), "/", i, sep = "")
    system(command)
  }
}

utils::globalVariables(c("palette_to_use", "full_palette",
                         "preview_info", "curated_dataset_names",
                         "dataset_selections", "default_celltype_curated",
                         "markdown_report_path", "all_zones"), "cytosel")

ggplot2::theme_set(cowplot::theme_cowplot())

options(shiny.maxRequestSize = 1000 * 200 * 1024 ^ 2, warn=-1,
        show.error.messages = FALSE)

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
#' @import forcats
#' @import ggplot2
#' @import sortable
#' @import scater
#' @import reactable
#' @import tidyverse
#' @import printr
#' @import promises
#' @importFrom rlang is_empty
#' @importFrom DT datatable
#' @importFrom clustifyr plot_gene
#' @importFrom readr write_lines read_tsv read_csv
#' @importFrom dplyr desc mutate_if distinct
#' @importFrom bsplus use_bs_popover shinyInput_label_embed shiny_iconlink bs_embed_tooltip use_bs_tooltip
#' @importFrom shinyjs useShinyjs hidden toggle reset
#' @importFrom grDevices dev.off pdf
#' @importFrom zip zip
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom SingleCellExperiment reducedDimNames reducedDims
#' @importFrom plotly plot_ly plotlyOutput renderPlotly layout hide_guides
#' @importFrom parallelly availableCores
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stringr str_split_fixed
#' @importFrom magrittr set_names
#' @importFrom shinydashboard box dashboardBody dashboardHeader dashboardSidebar
#' dashboardPage menuItem sidebarMenu sidebarMenuOutput tabItem tabItems
#' valueBoxOutput renderMenu updateTabItems tabBox
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom yaml read_yaml
#' @import rdrop2
#' @importFrom lubridate with_tz
#' @export
#' 
#' @param ... Additional arguments
cytosel <- function(...) {
  
  antibody_info <- cytosel_data$antibody_info |> dplyr::rename(Symbol = 
                                                  `Gene Name (Upper)`) |>
                    tidyr::drop_na()
  
  options(MulticoreParam=quote(MulticoreParam(workers=availableCores())))

  
  applications_parsed <- get_antibody_applications(antibody_info, 
                                                   'Symbol', 'Listed Applications')
  
  grch38 <- cytosel_data$grch38
  
  full_palette <- create_global_colour_palette()
  
  markdown_report_path <- system.file(file.path("report", "rmarkdown-report.Rmd"), 
                                      package = "cytosel")
  
  # devtools will find the file in the inst directory (move to top level)
  cytosel_token <- readRDS(system.file(file.path("token.rds"),
                                       package = "cytosel"))
  
  all_zones <- cytosel_data$time_zones
  
  ui <- tagList(
    includeCSS(system.file(file.path("www", "cytosel.css"),
                          package = "cytosel")),
    tags$head(
      tags$script(HTML("window.onbeforeunload = function() {return 'Your changes will be lost!';};"))
    ),
    tags$head(tags$style(".modal-dialog{ width:750px}")),
    tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),
    # styling the hover tooltips
    #  https://stackoverflow.com/questions/62360730/change-color-of-bstooltip-boxes-in-shiny
    tags$style(HTML("
                .tooltip > .tooltip-inner {
                width: 450px;
                color: white;
                background-color: black;
                }
                ")),
  
    # Use packages
    useShinyalert(force = TRUE),
    use_bs_tooltip(),
    useShinyjs(),
    
    dashboardPage(
    
    # Title
    dashboardHeader(title = "cytosel", tags$li(class = "dropdown", 
                                               htmlOutput("current_session_info",
                                    style = "width:100%; color:white; height:70%; margin-right:12px; margin-top:6.5px;
                                    margin-bottom:-2px;")),
                    tags$li(class = "dropdown", 
                                        actionLink("time_zone", "Set Time Zone",
                                        width = "90%", style = "margin-top:-0.5px;
                                        margin-right: 15px; margin-bottom:-2px;"))),
    
    dashboardSidebar(
      sidebarMenu(id = "tabs",
        # width = 300,
        menuItem("Get Started", tabName = "inputs", icon = icon("gear")),
        sidebarMenuOutput(outputId = 'output_menu'),
        menuItem("Documentation", tabName = "documentation", icon = icon("bookmark")),
        column(12, align = "center", offset = 0,
               hidden(div(id = "analysis_button", actionButton("start_analysis", "Run analysis!",
                           icon=icon("play", style = "color:black;"),
                           style = "margin-left: -5px;")))),
        column(12, align = "center", offset = 0,
               hidden(div(id = "download_button", downloadButton("downloadData",
               "Save panel", style = "color:black;
                           margin-top: 15px; margin-left:-10px; width: 60%;"))))
    )),
    dashboardBody(
      
      # https://stackoverflow.com/questions/52198452/how-to-change-the-background-color-of-the-shiny-dashboard-body
      tags$head(tags$style(HTML('
                                /* body */
                                .content-wrapper, .right-side {
                                background-color: #FFFFFF;
                                }
                                
                                #cell_cat_preview {
                                max-width: 30%;
                                }
                                
                                .fa {
                                 font-size: 15px;
                                }
                          
                                #base_inputs {
                                border: 2.5px dashed black;
                                }
                                
                                '))),
      
    # Tabs
    tabItems(
      tabItem("inputs",
              div(style = "margin-top: -10px;"),
              # fluidRow(column(3), tags$h3("Quickstart",
              #                             style = 'padding-left: 20px; padding-top: -27px;
              #                             pading-bottom: 20px')
              #          ),
              fluidRow(column(3, tags$h4(HTML("<em>Step 1: Select dataset.</em>"),
                                          style = 'padding-left: 16px; padding-top: 0px; 
                                          padding-bottom: 15px')),
                       column(1,
                              icon("circle-info", style = "margin-left: -11px; margin-top: 13px") %>%
                             bs_embed_tooltip(title = get_tooltip('input_scrnaseq'),
                                              placement = "right", html = "true")
              )),
              fluidRow(column(4, box(title = "Upload scRNA-seq data", status = "primary",
                                     width = 12, style = "padding-bottom: -15px",
                                     div(style = "margin-top: -30px; margin-bottom: -25px"),
        fileInput("input_scrnaseq",
                  label = "", accept = c(".rds")),
        div(style = "margin-top: -30px")
        )),
        tags$h4("or", style = "width = 80%; padding-left: -10px; padding-top: -40px"),
        div(style = "margin-top: -56px"),
        column(3, align = "center", style = "margin-left: 15px", box(title = "Select a curated dataset", status = "primary",
                      width = 12, actionButton("curated_dataset", "Browse datasets",
                                               style = "margin-top: -10px; margin-bottom: 10px;"
                                               ) 
                                     ))),
        fluidRow(column(3, tags$h4(HTML("<em> Step 2: Select a cell category to analyze.</em>"),
                          style = 'padding-left: 16px; padding-top: 
                          0px; padding-bottom: 15px;')),
                 column(1,
                        icon("circle-info", style = "margin-left: -10px; margin-top: 17px") %>%
                          bs_embed_tooltip(title = get_tooltip('coldata_column'),
                                           placement = "right", html = "true")
                 )
        ),
        fluidRow(column(4, box(title = "Browse dataset metadata",
                               div(style = "margin-top: -30px"),
                               selectInput("coldata_column", "", NULL, multiple=FALSE),
                               width = 12, status = "primary"))),
        textOutput("cell_cat_preview"),
        br(),
        # add padding space between elements
        fluidRow(column(4, bsCollapse(id = "advanced_collapse",
                       bsCollapsePanel(title = HTML(paste0(
                         "Advanced settings", tags$span(icon("sort-down",
                                                       style = "position:right; margin-left: 4px; margin-top: -4px;")))), style = "info",
                                    selectInput("assay",
                                     "Choose which assay to use",
                                     NULL),
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
                               bs_embed_tooltip(title = get_tooltip('applications'),
                               placement = "right", html = "true"))
                   ),
                   bsCollapsePanel(title = HTML(paste0(
                     "Upload previous analysis", tags$span(icon("sort-down",
                                                style = "position:right; margin-left: 4px; margin-top: -4px;")))), style = "info",
                     fileInput("read_back_analysis",
                               label = "Select a yml file from a previous run",
                               accept = c(".yml")) %>%
                       shinyInput_label_embed(
                         icon("circle-info") %>%
                           bs_embed_tooltip(title = get_tooltip('reupload'),
                                            placement = "right", html = "true")),
                     br()
                   )))),
        br(),
        br(),
        actionButton("create_reset", label = "Reset marker panel")
        ),
      
      tabItem("marker_selection",
                 # icon = icon("list"),
                 br(),
                 fluidRow(column(12, actionButton("markers_change_modal", "Add markers to panel"))),
                 hr(),
                 fluidRow(column(7, hidden(div(id = "marker_visualization",
                                               tags$span(icon("circle-info")
                                                         %>%
                                                           bs_embed_tooltip(title = get_tooltip('marker_visualization'),
                                                                            placement = "right", html = TRUE)))),),
                          column(5, hidden(div(id = "marker_display",
                   tags$span(icon("circle-info")
                           %>%
                   bs_embed_tooltip(title = get_tooltip('marker_display'),
                                    placement = "right", html = TRUE)))))),
                fluidRow(column(4, align = "center", div(style="display: inline-block; font-size: 15px", 
                                  htmlOutput("scratch_marker_counts"))),
                         column(4, align = "center", div(style="display: inline-block; font-size: 15px", 
                                  htmlOutput("selected_marker_counts")))),
                 fluidRow(column(8, uiOutput("BL"),
                                 style = "margin-bottom:0px;"
                                 ),
                          column(4, align = "left", plotOutput("legend", width = "400px",
                                               height = "500px"),
                                 style = " margin-top:-70px; margin-left: -50px;"))),
      tabItem("gene_expression",
              fluidRow(column(3,
                selectizeInput("genes_for_violin", "Select genes to view their expression:", NULL, multiple=T)),
                column(2, actionButton("add_violin_genes", "Add selected to panel", width = "100%"),
                       style = "margin-top:25px"),
                column(4, radioButtons("viol_viewer", label = "Select plot orientation",
                                       choices = c("By Marker", "By Cell Type"),
                                       selected="By Marker"))),
              br(),
              fluidRow(column(12, hidden(div(id = "viol_plot",
                    plotOutput("expression_violin", height="550px")))))
      ),
      tabItem("UMAP",
                 # icon = icon("globe"),
                 br(),
                 helpText(get_tooltip('umap')),
              br(),
              splitLayout(cellWidths = c(320, 280), selectInput("umap_options", 
                                "Select UMAP colouring:", choices = c("Cell Type",
                                                                    "Panel Marker"),
                                selected = "Cell Type", multiple=FALSE),
                              column(6, hidden(div(id = "umap_panel_cols",
                                         selectizeInput("umap_panel_options", 
                                                     "Color Heatmap by Panel Marker", 
                                                     choices = NULL, multiple = F,
                                                     width = "89%",
                                                     ))))),
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
             
                  splitLayout(cellWidths = c(250, 320),
                                          div(style="margin-right:10px; margin-bottom:25px;",
                              numericInput("n_genes", "Remove redundant genes", 
                                           value = 10, min = 1, width = "110%") %>%
                                shinyInput_label_embed(icon("circle-info") %>%
                              bs_embed_tooltip(title = get_tooltip('gene_removal'),
                                              placement = "right"))),
                          div(style="margin-left:10px; margin-top: 25px; margin-bottom:25px;",
                              actionButton("suggest_gene_removal", "View suggestions"))),
                 fluidRow(column(12, plotlyOutput("heatmap", height="550px"))),
                 br()
          ),

      tabItem("Metrics",
                 # icon = icon("ruler"),
                 br(),
                 helpText(get_tooltip('metrics')),
                 tags$span(icon("circle-info") %>%
                             bs_embed_tooltip(title = get_tooltip('metrics_explanation'),
                                              placement = "right", html = TRUE)),
                 # textOutput("cells_per_category"),
                fluidRow(
                column(8, plotlyOutput("metric_plot", height="550px"),),
                column(4, align = "center",
                       tabBox(width = NULL,
                  # Title can include an icon
                  tabPanel("Current Run Metrics",
                           DTOutput("current_run_metrics", width = "100%")
                  ),
                  tabPanel("Previous Run Metrics",
                           DTOutput("previous_run_metrics", width = "100%"))
                )
          ))),

      tabItem("alternative_markers",
                 # icon = icon("exchange-alt"),
                 br(),
                 helpText(get_tooltip("alternative_markers")),
                 splitLayout(cellWidths = c(180, 240),
                             div(selectInput("input_gene", "Input gene", 
                                             choices = NULL,
                                             selected = NULL,
                                             width = "100%")),
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
                fluidRow(column(12,
                                  br(),
                 reactableOutput("antibody_table"),
                # DT::dataTableOutput("antibody_table")))
          ))),
      tabItem("runs",
              tabsetPanel(type = "tabs",
                          tabPanel(uiOutput("current_run_name"), 
                                   fluidRow(column(6, h4("Run Parameters"),
                                                   DTOutput("summary_run_current")),
                                            column(6, h4("Run Metrics"),
                                                   DTOutput("metrics_run_current")))),
                          tabPanel(uiOutput("previous_run_1"),
                                   fluidRow(column(6, h4("Run Parameters"),
                                    DTOutput("summary_prev_1")),
                                      column(6, h4("Run Metrics"),
                                        DTOutput("metrics_run_prev_1")))),
                          tabPanel(uiOutput("previous_run_2"),
                                   fluidRow(column(6, h4("Run Parameters"),
                                                   DTOutput("summary_prev_2")),
                                            column(6, h4("Run Metrics"),
                                                   DTOutput("metrics_run_prev_2"))))
              )),
      tabItem("documentation",
              htmlOutput("cytosel_hyperlink")
              # h2("Documentation Preview:"),
              # htmlOutput("cytosel_preview"))
      )
    )
  )
  )
  )
  
  server <- function(input, output, session) {
    
    ### REACTIVE VARIABLES ###
    plots <- reactiveValues() # Save plots for download
    
    umap_top <- reactiveVal()
    umap_top_gene <- reactiveVal()
    umap_all_gene <- reactiveVal()
    umap_colouring <- reactiveVal()
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
    reupload_cell_types <- reactiveVal(FALSE)
    yaml <- reactiveVal()
    
    reset_panel <- reactiveVal(FALSE)
    valid_existing_panel <- reactiveVal(TRUE)
    default_category_curated <- reactiveVal()
    
    cell_types_missing_markers <- reactiveVal(NULL)
    time_zone_set <- reactiveVal()
    
    # run log variables #
    current_run_log <- reactiveVal()
    previous_run_log <- reactiveVal(NULL)
    previous_run_log_2 <- reactiveVal(NULL)
    
    proper_organism <- reactiveVal(TRUE)
    
    downloaded_content <- reactiveVal(FALSE)
    
    output$cytosel_hyperlink <-  renderUI({
      # url <- a("Cytosel Documentation", href="http://camlab-bioml.github.io/cytosel-doc/docs/intro")
      # tagList("Cytosel Documentation", url)
      tags$a(href="http://camlab-bioml.github.io/cytosel-doc/docs/intro",
             "Click Here to access the cytosel documentation in a new tab",
             target="_blank")
    })
    
    output$current_session_info <- renderUI({
      if (isTruthy(time_zone_set())) {
        time_zone <- time_zone_set()
      } else {
        time_zone <- Sys.timezone()
      }
      
      HTML(paste("session Id:", Sys.info()[["nodename"]],
            "<br/>", "server time zone:", time_zone))
    })
    
    observeEvent(input$time_zone, {
      output$current_time <- renderUI({
        HTML(paste("Approximate date and time in: ", input$time_zone_options,
                   "<br/>", "<b>",
                   as.character(with_tz(Sys.time(), 
                  tzone = input$time_zone_options)),
                  "<br>"))})
      showModal(time_zone_modal(all_zones, time_zone_set()))
    })
    
    observeEvent(input$pick_time_zone, {
      req(input$time_zone_options)
      
      Sys.setenv(TZ=input$time_zone_options)
      time_zone_set(input$time_zone_options)
      removeModal()
    })
    
    
    # output$cytosel_preview <- renderUI({
    #   tags$iframe(src="http://camlab-bioml.github.io/cytosel-doc/docs/intro", height=600, width=1100)
    # })
    
    observeEvent(input$curated_dataset, {
      
      showModal(curated_dataset_modal(curated_dataset_names))
    })
    
    observeEvent(input$curated_options, {
    
      output$curated_set_preview <- renderPrint({HTML(preview_info[names(preview_info) ==
                                                              input$curated_options])})
    })
    
    observeEvent(input$pick_curated, {
      req(input$curated_options)
      
      future_promise({
      removeModal()
      pre_upload_configuration()
      
      withProgress(message = 'Configuring curated selection', value = 0, {
        setProgress(value = 0)
      
      if (!file.exists(file.path(tempdir(), dataset_selections[names(dataset_selections) ==
                                                        input$curated_options]))) {
          incProgress(detail = "Downloading curated dataset")
          rdrop2::drop_download(paste("cytosel/", dataset_selections[names(dataset_selections) ==
                                                                       input$curated_options],
                                sep = ""),
                                local_path = tempdir(),
                                overwrite = T,
                                dtoken = cytosel_token)
      }
        
        setProgress(value = 0.5)
        incProgress(detail = "Reading input dataset")
        
        input_sce <- read_input_scrnaseq(file.path(tempdir(), 
                                          dataset_selections[names(dataset_selections) ==
                                          input$curated_options]))
        
        incProgress(detail = "Parsing gene names and assays")
        
        setProgress(value = 1)
      })
      
      default_category_curated(default_celltype_curated[names(default_celltype_curated) ==
                           input$curated_options])
      
      post_upload_configuration(input_sce)
    })
    })
    
    
    ### UPLOAD FILE ###
    observeEvent(input$input_scrnaseq, {
      
      future_promise({
      default_category_curated(NULL)
      withProgress(message = 'Configuring input selection', value = 0, {
      setProgress(value = 0)
      pre_upload_configuration()
      setProgress(value = 0.25)
      incProgress(detail = "Reading input dataset")
      input_sce <- read_input_scrnaseq(input$input_scrnaseq$datapath)
      incProgress(detail = "Parsing gene names and assays")
      setProgress(value = 0.85)
      found_human <- check_for_human_genes(input_sce)
      if (!isTruthy(found_human)) {
        throw_error_or_warning(type = "error", message = "cytosel has detected that the gene names 
                               provided are not human. Currently, only human datasets are supported.")
        proper_organism(found_human)
      }
      setProgress(value = 1)
      })
      
      req(proper_organism())
      post_upload_configuration(input_sce)
    })
    })
    
    observeEvent(input$coldata_column, {
      req(sce())
      column(input$coldata_column)
      
      len_possible_cats <- length(unique(sce()[[input$coldata_column]]))
      num_limit <- ifelse(len_possible_cats <= 3, len_possible_cats, 3)
      cats_to_show <- unique(sce()[[input$coldata_column]])[1:num_limit]
      others_addition <- ifelse(len_possible_cats <= 3, "", "and others")
      var_group <- ifelse(len_possible_cats < 2, "grouping", "groupings")
      
      output$cell_cat_preview <- renderText({paste(len_possible_cats,
                                                   var_group, "in selected category, including",
                                                   paste(cats_to_show,
                                                         collapse = ", "),
                                                   others_addition,
                                                   sep = " ")})
      
      if (!isTruthy(reupload_cell_types())) {
        updateSelectInput(session, "user_selected_cells",
                          unique(sce()[[input$coldata_column]]))
        specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      }
      
      cell_types_to_keep(NULL)

      reupload_analysis(FALSE)
      reupload_cell_types(FALSE)
      
      
      })
    
    observeEvent(input$add_selected_to_analysis, {
      req(sce())
      req(input$coldata_column)
      specific_cell_types_selected(input$user_selected_cells)
      
      removeModal()
      
      if (length(specific_cell_types_selected()) > 0) {
        showNotification("Setting subset to select cell types",
                         duration = 3)
      } else {
        showNotification("Empty subset selection, defaulting to using all cell types in analysis.",
                         duration = 3)
        specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      }
      
    })
    
    observeEvent(input$user_selected_cells, {
      req(sce())
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
    
    observeEvent(input$assay, {
      pref_assay(input$assay)
    })
    
    ### bring up modal to add markers ###
    observeEvent(input$markers_change_modal, {
      req(sce())
      req(current_markers())
      req(allowed_genes())
      req(fms())
    
      showModal(markers_add_modal(allowed_genes(), names(fms()[[1]])))
    })
    
    # ### ANTIBODY EXPLORER ###
    output$antibody_table <- renderReactable({
      req(current_markers())
      req(df_antibody())

      reactable(df_antibody(),
                searchable = TRUE,
                filterable = TRUE,
                columns = list(`Host Species` = colDef(
                  filterInput = function(values, name) {
                    tags$select(
                      # Set to undefined to clear the filter
                      onchange = sprintf("Reactable.setFilter('antibody-select', '%s', event.target.value || undefined)", name),
                      # "All" has an empty value to clear the filter, and is the default option
                      tags$option(value = "", "All"),
                      lapply(unique(values), tags$option),
                      "aria-label" = sprintf("Filter %s", name),
                      style = "width: 100%; height: 28px;"
                    )
                  }
                ),
                `Product Category Tier 3` = colDef(
                  filterInput = function(values, name) {
                    tags$select(
                      # Set to undefined to clear the filter
                      onchange = sprintf("Reactable.setFilter('antibody-select', '%s', event.target.value || undefined)", name),
                      # "All" has an empty value to clear the filter, and is the default option
                      tags$option(value = "", "All"),
                      lapply(unique(values), tags$option),
                      "aria-label" = sprintf("Filter %s", name),
                      style = "width: 100%; height: 28px;"
                    )
                  }
                ),
                `External Link` = colDef(html = T)
                ),
                sortable = TRUE,
                elementId = "antibody-select")
    })

    ### PLOTS ###
    output$all_plot <- renderPlotly({
      req(umap_all())
      req(column())
      req(plots$all_plot)
      
      columns <- column()

      plots$all_plot
    })
    
    output$top_plot <- renderPlotly({
      req(umap_top())
      req(column())
      req(plots$top_plot)
      
      plots$top_plot
      
    })
    
    output$heatmap <- renderPlotly({
      req(heatmap())
      req(column())
      
     heatmap()
    })
    
    output$metric_plot <- renderPlotly({
      req(metrics())
      req(cells_per_type())
      req(plots$metric_plot)
     
      plots$metric_plot
      
    })
    
    output$current_run_metrics <- renderDT({
      req(current_metrics())
      current_metrics()$summary
      })

    output$previous_run_metrics <- renderDT({
        req(previous_metrics())
        previous_metrics()$summary})
    
    
    ##### Current run logs #####
    
    output$current_run_name <- renderText({
      req(current_markers())
      req(current_run_log())
      as.character(current_run_log()$map$`Time`)
    })


    output$summary_run_current <- renderDT({
      req(current_metrics())
      req(current_run_log())
      current_run_log()$frame}, server = TRUE)

    output$metrics_run_current <- renderDT({
      req(current_metrics())
      req(current_run_log())
      current_run_log()$metrics}, server = TRUE)
    
    ##### Previous run logs #####

    output$previous_run_1 <- renderText({
        req(previous_run_log())
        as.character(previous_run_log()$map$`Time`)
      })

    output$summary_prev_1 <- renderDT({
      req(previous_run_log())
      previous_run_log()$frame},  server = TRUE)

    output$metrics_run_prev_1 <- renderDT({
      req(previous_run_log())
      previous_run_log()$metrics}, server = TRUE)
    
    output$previous_run_2 <- renderText({
      req(previous_run_log_2())
      as.character(previous_run_log_2()$map$`Time`)
    })
    
    output$summary_prev_2 <- renderDT({
      req(previous_run_log_2())
      previous_run_log_2()$frame},  server = TRUE)
    
    output$metrics_run_prev_2 <- renderDT({
      req(previous_run_log_2())
      previous_run_log_2()$metrics}, server = TRUE)
    
    ### ANALYSIS ###
    observeEvent(input$start_analysis, {
      req(input$panel_size)
      req(input$coldata_column)
      req(sce())
      
      if (!is_empty(input$bl_top)) {
        current_panel_not_valid <- setdiff(input$bl_top, rownames(sce()))
        if (length(current_panel_not_valid) > 0) {
          valid_existing_panel(FALSE)
          current_pan_not_valid_modal(current_panel_not_valid)
        }
      }
    
      if (!isTruthy(specific_cell_types_selected()) | length(specific_cell_types_selected()) < 1) {
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
      
        req(valid_existing_panel())
        req(proceed_with_analysis())
        
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
      
      
      withProgress(message = 'Configuring analysis', value = 0, {
        incProgress(detail = "Acquiring data")
        req(proceed_with_analysis())
        req(any_cells_present())
        
          ## Set initial markers:
          scratch_markers_to_keep <- input$bl_scratch
          
          columns <- good_col(sce()[,sce()$keep_for_analysis == "Yes"], input$coldata_column)
          column(columns$good)
          col <- columns$bad
          
          
          if (input$subsample_sce) {
            incProgress(detail = "Subsampling data")
            
            avail_cols <- which(sce()$keep_for_analysis == "Yes")
            to_subsample <- sample(avail_cols, min(length(avail_cols), 2000),
                                   replace = FALSE)
            
            sce(create_keep_vector_during_subsetting(sce(), to_subsample))
            
            cell_dist_subsample <- create_table_of_hetero_cat(sce()[,sce()$keep_for_analysis == "Yes"],
                                                              input$coldata_column) |>
              # filter(Freq > 0 & (Freq < 2 | `Proportion Percentage` < 0.5))
              filter(Freq > 0 & Freq < 2)
            
            if (nrow(cell_dist_subsample) > 0) {
              subsampling_error_modal(unique(cell_dist_subsample[,1]))
              proceed_with_analysis(FALSE)
            }
            
          } 
          
          setProgress(value = 1)
          
          withProgress(message = 'Starting computations', value = 0, {
            req(proceed_with_analysis())
            
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
            
            setProgress(value = 0.5)
            
            incProgress(detail = "Computing cell type markers")
            
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
            
            setProgress(value = 0.75)
            
            if(!is_empty(input$bl_top)) {
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
                markers_list <- get_markers(fms(), 
                                       # Adding 10 to make sure panel size is approximate 
                                       # since a) the same marker is selected multiple times and
                                       # b) excess markers are removed
                                       input$panel_size, 
                                       input$marker_strategy, 
                                       sce()[,sce()$keep_for_analysis == "Yes"],
                                       allowed_genes())
                
                if (isTruthy(markers_list$missing)) {
                  throw_error_or_warning(type = 'notification',
                                         message = paste("No markers were found for the following cell types: ",
                                                         paste(markers_list$missing, 
                                                               collapse = ", "), 
                                                         ". This is likely because there are too few cells of these types."),
                                         duration = 10)
                  cell_types_missing_markers(markers_list$missing)
                  
                }
                
                markers <- markers_list$marker[!is.na(markers_list$marker)]
                
                setProgress(value = 1)
                
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
                                        current_markers()$top_markers) |>
                          mutate(`Host Species` = factor(`Host Species`),
                                 `Product Category Tier 3` = factor(`Product Category Tier 3`),
                                 `KO Status` = factor(`KO Status`),
                                 `Clone Number` = factor(`Clone Number`),
                                 `External Link` = paste0('<a href="',`Datasheet URL`, '"',
                                                          ' target="_blank" rel="noopener noreferrer"',
                                                          '>',"View in Abcam website",'</a>')) |>
                          dplyr::select(-c(`Datasheet URL`)))
            
            update_analysis()
            
            updateSelectizeInput(
              session = session,
              inputId = "umap_panel_options",
              # choices = current_markers()$top_markers,
              # selected = current_markers()$top_markers[1])
              choices = allowed_genes(),
              selected = current_markers()$top_markers[1])
            
            updateSelectizeInput(
              session = session,
              inputId = "genes_for_violin",
              choices = allowed_genes(),
              server = T)

            if (!isTruthy(first_render_outputs()) & isTruthy(current_markers())) {
              
              output$output_menu <- renderMenu(expr = {
                sidebarMenu(
                  menuItem("Marker selection", tabName = "marker_selection", 
                           icon = icon("barcode")),
                  menuItem("Gene Expression", tabName = "gene_expression", 
                           icon = icon("arrows-alt")),
                  menuItem("UMAP", tabName = "UMAP", 
                           icon = icon("circle-nodes")),
                  menuItem("Heatmap", tabName = "Heatmap", 
                           icon = icon("th-large")),
                  menuItem("Metrics", tabName = "Metrics", 
                           icon = icon("line-chart")),
                  menuItem("Alternative Markers", tabName = "alternative_markers", 
                           icon = icon("arrows-h")),
                  menuItem("Antibody Explorer", tabName = "antibody_explorer", 
                           icon = icon("list-alt")),
                  menuItem("Run Logs", tabName = "runs", 
                           icon = icon("magnifying-glass"))
                )
              })
              
              first_render_outputs(TRUE)
              setProgress(value = 1)
            }
            reupload_analysis(FALSE)
            updateTabsetPanel(session, "tabs", "marker_selection")
            
          } else {
            unique_element_modal(col)
          }
          
          
          
          })
          
      })
      
      # reset downloaded status when analysis is rerun
      downloaded_content(FALSE)
      toggle(id = "download_button", condition = isTruthy(current_markers()))
      
    })
    
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
      
      output$selected_marker_counts <- renderText({paste("<B>",
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
    
    
    #### Add genes from expression tab #####
    
    observeEvent(input$add_violin_genes, {
      req(input$genes_for_violin)
      
      if(!is.null(input$genes_for_violin) && 
         all(input$genes_for_violin %in% allowed_genes()) &&
         !all(input$genes_for_violin %in% input$bl_scratch) &&
         !all(input$genes_for_violin %in% input$bl_top)) {
        
        cm <- current_markers()
        markers <- list(recommended_markers = cm$recommended_markers,
                        scratch_markers = input$bl_scratch,
                        top_markers = unique(c(input$genes_for_violin,
                        setdiff(cm$top_markers, input$bl_scratch))))
        
        # SMH
        current_markers(
          set_current_markers_safely(markers, fms())
        )
        
        num_markers_in_selected(length(current_markers()$top_markers))
        num_markers_in_scratch(length(current_markers()$scratch_markers))
        
        update_BL(current_markers(), num_markers_in_selected(),
                  num_markers_in_scratch(),
                  names(fms()[[1]]))
        
      } 
      updateTabsetPanel(session, "tabs", "marker_selection")
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
      
      suggestions(suggest_genes_to_remove(cmat, input$n_genes))
      
      showModal(suggestion_modal(suggestions = suggestions(),
                                 possible_removal = current_markers()$top_markers))
    })
    
    observeEvent(input$remove_suggested, { # Remove suggested markers
      req(current_markers())
      cm <- current_markers()
      
      if (!is.null(input$markers_to_remove)) {
        remove_marker <- c(input$markers_to_remove)
        # num_markers_in_selected(length(current_markers()$top_markers))
        # num_markers_in_scratch(length(current_markers()$scratch_markers))

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
        
        withProgress(message = 'Updating panel', value = 0, {
          incProgress(6, detail = "Computing alternatives")
          
          # Make this sampling dependent on the input sample argument
          # exclude genes that are already in the marker panel
          replacements(
            compute_alternatives(input$input_gene, 
                                 sce()[,sce()$keep_for_analysis == "Yes"], 
                                 pref_assay(), input$number_correlations,
                                 allowed_genes()[!allowed_genes() %in% current_markers()$top_markers &
                                                   !allowed_genes() %in% current_markers()$scratch_markers]) |>
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
      
    observeEvent(input$send_markers, {# Send alternative markers to selected markers panel
    
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
      
      yaml_back <- read_back_in_saved_yaml(input$read_back_analysis$datapath)
      
      if (isTruthy(yaml_back$`Target panel size`)) {
        updateNumericInput(session, "panel_size", value = yaml_back$`Target panel size`)
      }
      
      markers_reupload(list(top_markers = yaml_back$`Selected marker panel`,
      scratch_markers = yaml_back$`Scratch marker panel`))
      if (is_empty(input$bl_top)) {
        showNotification("Updating the current panel to the markers from the reuploaded yml.",
                         type = "message", duration = 4)
      }
      
      yaml_length <- sum(c(length(markers_reupload()$top_markers),
                           length(markers_reupload()$scratch_markers)))
      if (yaml_length != yaml_back$`Target panel size` & yaml_length > 0) {
        updateNumericInput(session, "panel_size", value = yaml_length)
      }
      
      if (isTruthy(sce())) {
        if (yaml_back$`Number of columns (cells)` != ncol(sce()) |
            yaml_back$`Number of rows (features)` != nrow(sce())) {
          reupload_warning_modal("Incompatible dimensions", "Number of cells and genes")
    }
        if (yaml_back$`Heterogeneity source` %in% colnames(colData(sce()))) {
        updateSelectInput(session, "coldata_column", choices = colnames(colData(sce())),
        selected = yaml_back$`Heterogeneity source`)
        updateSelectInput(session, "user_selected_cells", 
                          yaml_back$`Cell Types Analyzed`)
        specific_cell_types_selected(yaml_back$`Cell Types Analyzed`)
        # cell_types_to_keep(yaml_back$`User selected cells`)
        reupload_cell_types(TRUE)
        
        } else {
          reupload_warning_modal("Cell type not found", yaml_back$`Heterogeneity source`)
        }
        updateCheckboxInput(session, inputId = "precomputed_dim",
                            value = yaml_back$`Pre-computed UMAP`)
        use_precomputed_umap(yaml_back$`Pre-computed UMAP`)
        
        if (isTruthy(use_precomputed_umap())) {
          possible_umap_dims(detect_umap_dims_in_sce(sce()))
          showModal(pre_computed_umap_modal(possible_umap_dims()))
        }
      }
      
    if (isTruthy(yaml_back$`Assay used`)) {
      pref_assay(yaml_back$`Assay used`)
    }
    
    if (isTruthy(yaml_back$`Min Cell Category cutoff`)) {
      cell_min_threshold(yaml_back$`Min Cell Category cutoff`)
      updateNumericInput(session, "min_category_count", value = yaml_back$`Min Cell Category cutoff`)
    }
    
    updateCheckboxInput(session, "subsample_sce", value = yaml_back$`Subsampling Used`)
    
    if (isTruthy(yaml_back$`Marker strategy`)) {
      updateRadioButtons(session, "marker_strategy", selected = yaml_back$`Marker strategy`)
    }
    
    if (isTruthy(yaml_back$`Marker strategy`)) {
      updateRadioButtons(session, "marker_strategy", selected = yaml_back$`Marker strategy`)
    }
    
    
    if (isTruthy(yaml_back$`Antibody applications`)) {
      updateSelectInput(session, "select_aa", selected = yaml_back$`Antibody applications`)
    }
    
    reupload_analysis(TRUE)
    
    })
    
    
    observeEvent(input$create_reset, {
      req(sce())
      req(current_markers())
      
      showModal(reset_analysis_modal())
    })
    
    observeEvent(input$reset_marker_panel, {
      
      reset_panel(TRUE)
      current_markers(list(top_markers = NULL, recommended_markers = NULL,
                           scratch_markers = NULL))
      reupload_analysis(FALSE)
      markers_reupload(NULL)
      num_markers_in_selected(0)
      num_markers_in_scratch(0)
      reset("bl_top")
      reset("bl_scratch")
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                names(fms()[[1]]))
      removeModal()
      updateTabsetPanel(session, "tabs", "marker_selection")
      valid_existing_panel(TRUE)
    
      })
    
    observeEvent(input$reset_marker_panel_reupload, {
      showModal(reset_analysis_modal())
      
    })
    
    observeEvent(input$umap_options, {
      req(sce())
      
      umap_colouring(input$umap_options)
      
      toggle(id = "umap_panel_cols", condition = input$umap_options != "Cell Type")
      
      observeEvent(input$umap_panel_options, {
        req(input$umap_panel_options)
        
        if (umap_colouring() == "Cell Type") {
          plots$all_plot <- suppressWarnings(plot_ly(umap_all(), x=~UMAP_1, y=~UMAP_2, color=~get(input$coldata_column),
                                                     text=~get(input$coldata_column), 
                                                     type='scatter', hoverinfo="text", colors=cytosel_palette()) %>% 
                                               layout(title = "UMAP all genes", showlegend = F))
          
          plots$top_plot <- suppressWarnings(plot_ly(umap_top(), x=~UMAP_1, y=~UMAP_2, color=~get(input$coldata_column),
                                                     text=~get(input$coldata_column),
                                                     type='scatter', hoverinfo="text", colors=cytosel_palette()) %>%
                                               layout(title = "UMAP selected markers"))
          umap_top_gene(FALSE)
          umap_all_gene(FALSE)
          
        } else {
          
          gene_plot_top <- plot_gene(assay(sce()[,sce()$keep_for_analysis == "Yes"], pref_assay()),
                                     umap_top() |> select(UMAP_1, UMAP_2) |> drop_na(),
                                     input$umap_panel_options)
          
          gene_plot_all <- plot_gene(assay(sce()[,sce()$keep_for_analysis == "Yes"], pref_assay()),
                                     umap_all() |> select(UMAP_1, UMAP_2) |> drop_na(),
                                     input$umap_panel_options)
          
          
          umap_top_gene(as.data.frame(gene_plot_top[[1]]$data))
          umap_all_gene(as.data.frame(gene_plot_all[[1]]$data))
          
          umap_top_gene(merge(umap_top_gene(), umap_top() |> select(input$coldata_column),
                              by.x = "cell", by.y = "row.names") |> 
                          mutate(lab = paste(get(input$coldata_column), ": ",
                                             "\n",
                                             round(get(input$umap_panel_options), 3), sep = "")) |>
                          rename(Expression = input$umap_panel_options))
          
          umap_all_gene(merge(umap_all_gene(), umap_all() |> select(input$coldata_column),
                              by.x = "cell", by.y = "row.names") |> 
                          mutate(lab = paste(get(input$coldata_column), ": ",
                                             "\n",
                                             round(get(input$umap_panel_options), 3), sep = "")) |>
                          rename(Expression = input$umap_panel_options))
          
          plots$top_plot <- suppressWarnings(plot_ly(umap_top_gene(),
                                    x = ~UMAP_1, y = ~UMAP_2, 
                                    color = ~Expression,
                                    type='scatter',
                                    text = ~lab,
                                    hoverinfo = "text",
                                    colors = c("grey60",
                                               "blue")) %>%
            # cytosel_palette()[current_markers()$associated_cell_types[input$umap_panel_options]])) %>%
            layout(title = "UMAP selected markers"))
          
          plots$all_plot <- hide_guides(suppressWarnings(plot_ly(umap_all_gene(),
                                    x = ~UMAP_1, y = ~UMAP_2, 
                                    color = ~Expression,
                                    type='scatter',
                                    text = ~lab,
                                    hoverinfo = "text",
                                    colors = c("grey60",
                                               "blue")) %>%
                                                 # cytosel_palette()[current_markers()$associated_cell_types[input$umap_panel_options]])) %>%
            layout(title = "UMAP all genes",
                   showlegend = F)))
          
        }
      })
      
    })
    
    to_listen_violin <- reactive({
      list(input$genes_for_violin, input$viol_viewer)
    })
    
    observeEvent(to_listen_violin(), {
      req(sce())
      
      observe(toggle(id = "viol_plot", condition = isTruthy(input$genes_for_violin)))

      if (isTruthy(input$genes_for_violin)) {
        
        viol_frame <- suppressWarnings(
          make_violin_plot(sce()[,sce()$keep_for_analysis == "Yes"],
                           input$genes_for_violin, input$coldata_column, pref_assay()))
        
        if (input$viol_viewer == "By Marker") {
          output$expression_violin <- renderPlot({
            suppressWarnings(ggplot(viol_frame, aes(x = `Cell Type`, y = Expression, 
                                                    fill = `Cell Type`)) + geom_violin(alpha = 1) +
                               theme(axis.text.x = element_text(angle = 90)) +
                               facet_wrap(~Gene, scales="free_y") + 
                               scale_fill_manual(values = cytosel_palette()))
          })
        } else {
          output$expression_violin <- renderPlot({
            suppressWarnings(ggplot(viol_frame, aes(x = `Cell Type`, y = Expression, 
                                                    fill = Gene)) + geom_violin(alpha = 1) +
                               theme(axis.text.x = element_text(angle = 90)) +
                               # facet_wrap(~Gene) + 
                               scale_fill_manual(values = full_palette))
          })
          
        }
        
      } else {
        output$expression_violin <- renderPlot({
        ggplot() + theme_void()
        })
      }
    }, ignoreNULL = FALSE)
    
    
    ### UPDATE SELECTED MARKERS ###
    update_BL <- function(markers, top_size, scratch_size, unique_cell_types,
                          adding_alternative = F) {
      set.seed(12345L)
      
      unique_cell_types <- sort(unique_cell_types)
      
      palette_to_use <- full_palette[1:length(unique_cell_types)]
      names(palette_to_use) <- unique_cell_types
      
      cytosel_palette(palette_to_use)
      
      if (length(markers$top_markers) > 0) {
        markers$top_markers <- sort(markers$top_markers)
      }
      
      if (length(markers$scratch_markers) > 0) {
        markers$scratch_markers <- sort(markers$scratch_markers)
      }
      
      labels_top <- lapply(markers$top_markers, 
                           function(x) div(x, map_gene_name_to_antibody_icon(x, markers), style=paste('padding: 3px; color:', 
                                                                                                      set_text_colour_based_on_background(cytosel_palette()[ markers$associated_cell_types[x]]), '; background-color:', 
                                                                                                   cytosel_palette()[ markers$associated_cell_types[x] ])))
      labels_scratch <- lapply(markers$scratch_markers, 
                               function(x) div(x, map_gene_name_to_antibody_icon(x, markers), style=paste('padding: 3px; color:', 
                                                                                                          set_text_colour_based_on_background(cytosel_palette()[ markers$associated_cell_types[x]]), '; background-color:', 
     
                                                                                   cytosel_palette()[ markers$associated_cell_types[x] ])))
 
      # output$legend <- renderPlot(cowplot::ggdraw(get_legend(cytosel_palette())))
      
      output$legend <- renderPlot({
        op <- par(family = "sans", mar = c(3, 2, 3, 3))
        plot(NULL ,xaxt='n',yaxt='n',bty='n',
                                       ylab='',xlab='', xlim=0:1, ylim=0:1)
      legend("top", legend = names(cytosel_palette()), 
             pch=15, pt.cex=2.2, cex=1.2, bty='n',
             ncol = ceiling(length(cytosel_palette())/15),
             col = cytosel_palette())
      mtext("Cell Type", cex=1.4)
      par(op)})
     
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
    
    #### Nested app functions ####
    
    ### UPDATE ANALYSIS ###
    update_analysis <- function() {
      
      withProgress(message = 'Updating visualizations', value = 0, {

        
        ## Re-set the set of allowed genes (these may have changed if a different
        ## antibody application is selected)
        allowed_genes(
          get_allowed_genes(input$select_aa, applications_parsed, 
                            sce()[,sce()$keep_for_analysis == "Yes"])
        )
        
        updateSelectInput(session = session,
        inputId = "input_gene",
        choices = current_markers()$top_markers,
        selected = NULL)
        
        num_markers_in_selected(length(current_markers()$top_markers))
        num_markers_in_scratch(length(current_markers()$scratch_markers))
        
        update_BL(current_markers(), num_markers_in_selected(),
                  num_markers_in_scratch(),
                  names(fms()[[1]]))
        
        setProgress(value = 0.25)
        
        
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
        
        plots$all_plot <- suppressWarnings(plot_ly(umap_all(), x=~UMAP_1, y=~UMAP_2, color=~get(columns[1]), text=~get(columns[1]), 
                                 type='scatter', hoverinfo="text", colors=cytosel_palette()) %>% 
          layout(title = "UMAP all genes", showlegend = F))
        
        plots$top_plot <- suppressWarnings(plot_ly(umap_top(), x=~UMAP_1, y=~UMAP_2, color=~get(columns[1]), text=~get(columns[1]), 
                                 type='scatter', hoverinfo="text", colors=cytosel_palette()) %>% 
          layout(title = "UMAP selected markers"))
        
        setProgress(value = 0.5)
        
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
        
        setProgress(value = 0.75)
        
        # Update metrics
        incProgress(detail = "Computing panel score")
        scores <- get_scores(sce()[,sce()$keep_for_analysis == "Yes"], 
                             column(), current_markers()$top_markers, pref_assay())
        metrics(scores)
        
        # update previous metrics before current
        
        if (isTruthy(current_metrics())) {
          previous_metrics(list(all = current_metrics()$all |> mutate(Run = "Previous Run"),
                                summary = current_metrics()$summary))
        }
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
        m$Counts <- plyr::mapvalues(as.character(m$what),
                                    from = names(cpt), 
                                    to = cpt)
        # m$what_tally <- plyr::mapvalues(as.character(m$what),
        #                           from = names(cpt), 
        #                           to = paste0(names(cpt), " (n = ", cpt, ")"))
        # m$what_tally <- as.factor(m$what_tally)
        # 
        # m$what_tally <- fct_reorder(m$what_tally, desc(m$score))
        # m$what_tally <- fct_relevel(m$what_tally, 
        #                             paste0("Overall (n = ", cpt['Overall'], ")"))
        # m$what_tally <- fct_rev(m$what_tally)
        
        # m[is.na(m)] <- 0
        
        m <- m |> drop_na(score)
        
        current_metrics(m |> mutate(Run = "Current Run"))
        current_metrics(list(all = current_metrics(),
                             summary = current_metrics() |>
                               mutate(Counts = as.numeric(Counts)) |>
                               group_by(what, Counts) |>
                               summarize(mean_score = round(mean(score),
                                                            3)) |>
                               arrange(desc(Counts)) |> rename(`Cell Type` = what,
                                                               `Mean Score` = mean_score)))
        
        if (isTruthy(previous_metrics()$all)) {
          all_metrics <- rbind(previous_metrics()$all, current_metrics()$all) |> 
            mutate(what = factor(what,
                                  levels = c(rev(unique(what[what != "Overall"])), "Overall")))
        } else {
          all_metrics <- current_metrics()$all |>  mutate(what = factor(what,
                    levels = c(rev(unique(what[what != "Overall"])), "Overall")))
        }
        
        plots$metric_plot <- suppressWarnings(plot_ly(all_metrics, x = ~score, y = ~what, 
                                     color = ~Run,
                                     type='box', hoverinfo = 'none') %>% 
          layout(boxmode = "group",
                 xaxis = list(title="Score"),
                 yaxis = list(title="Source")))
        
        previous_run_log_2(previous_run_log())
        previous_run_log(current_run_log())

        run_config <- create_run_param_list(marker_list = current_markers(),
                                            input_file = input$input_scrnaseq$datapath,
                                            assay_used = pref_assay(),
                                            het_source = column(),
                                            panel_size = input$panel_size,
                                            cell_cutoff_value = as.integer(cell_min_threshold()),
                                            subsample = input$subsample_sce,
                                            marker_strat = input$marker_strategy,
                                            antibody_apps = input$select_aa,
                                            selected_cell_types = cell_types_to_keep(),
                                            precomputed_umap_used = input$precomputed_dim,
                                            num_cells = ncol(sce()),
                                            num_genes = nrow(sce()),
                                            metrics = current_metrics()$summary)

        config_df <- tibble::enframe(run_config) %>%
          dplyr::mutate(value = purrr::map_chr(value, toString)) |>
          `colnames<-`(c("Parameter", "Value")) |>
          filter(! Parameter %in% c("Run Metrics", "Input file"))
        
        current_run_log(list(map = run_config,
                             frame = config_df,
                             metrics = current_metrics()$summary))
        
        # Show help text popover
        shinyjs::show(id = "marker_visualization")
        shinyjs::show(id = "marker_display")
        
        setProgress(value = 1)
        
        
      })
    }
    
    pre_upload_configuration <- function() {
      updateCheckboxInput(inputId = "precomputed_dim", value = F)
      use_precomputed_umap(FALSE)
      umap_precomputed_col(NULL)
    }
    
    post_upload_configuration <- function(input_sce) {
      input_sce <- detect_assay_and_create_logcounts(input_sce)
      input_sce <- parse_gene_names(input_sce, grch38)
      sce(input_sce)
      
      shinyjs::hide(id = "download_button")
      
      toggle(id = "analysis_button", condition = isTruthy(sce()))
      
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
          showNotification("Setting assay to logcounts. This can be changed under Advanced settings",
                           type = 'message',
                           duration = 4)
        } else {
          showNotification(paste("Setting assay to ", input_assays[1],
                                 ". This can be changed under Advanced settings",sep = ""),
                           type = 'message',
                           duration = 4)
        }
      }
      
      updateSelectInput(
        session = session,
        inputId = "assay",
        choices = input_assays,
        selected = input_assays[1]
      )
      
      selection <- ifelse(isTruthy(default_category_curated()),
                          default_category_curated(),
                          colnames(colData(sce()))[!grepl("keep_for_analysis", 
                                                          colnames(colData(sce())))][1])
      
      
      updateSelectInput(
        session = session,
        inputId = "coldata_column",
        choices = colnames(colData(sce()))[!grepl("keep_for_analysis", 
                                                  colnames(colData(sce())))],
        selected = selection
      )
      
      if (!isTruthy(input$coldata_column)) {
        column(colnames(colData(sce()))[1])
      } else {
        column(input$coldata_column)
      }
      
      possible_umap_dims(detect_umap_dims_in_sce(sce()))
      
      toggle(id = "precomputed", condition = length(possible_umap_dims()) > 0)
      
      if (isTruthy(input$bl_top) | isTruthy(current_markers())) {
        showModal(reset_option_on_upload_modal())
      }
    }
    
    
    ### SAVE PANEL ###
    output$downloadData <- downloadHandler(
      filename <- paste0("Cytosel-Panel-", Sys.Date(), ".zip"),
  
      content = function(fname) {
        showNotification("Rendering output report and config file, this may take a few moments..",
                         duration = 7)
        future_promise({download_data(fname,
                      current_run_log()$map, plots, heatmap(), 
                      df_antibody(), markdown_report_path, current_metrics()$summary)})
      },
      contentType = "application/zip"
    )
    
  }
  
  shinyApp(ui, server, ...)
}


