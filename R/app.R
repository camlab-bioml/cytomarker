
curated_datasets <- utils::read.delim(system.file("ts_datasets.tsv", package = "cytomarker"),
                                      sep = "\t")

other_curated <- utils::read.delim(system.file("other_datasets.tsv", package = "cytomarker"),
                                   sep = "\t")

for (i in curated_datasets$tissue) {
  if (file.exists(file.path(tempdir(), "/", paste(i, ".rds", sep = "")))) {
    command <- paste('rm ', tempdir(), "/", paste(i, ".rds", sep = ""), sep = "")
    system(command)
  }
}

utils::globalVariables(c("palette_to_use", "full_palette",
                         "preview_info", "curated_datasets", "compartments",
                         "dataset_labels", "other_dataset_labels", 
                         "all_curated_dataset_labels",
                         "antibody_info", "markdown_report_path", 
                         "all_zones"), "cytomarker")

options(shiny.maxRequestSize = 1000 * 200 * 1024 ^ 2, warn=-1,
        show.error.messages = FALSE)

yaml <- read_yaml(system.file("config.yml", package = "cytomarker"))
USE_ANALYTICS <- yaml$use_google_analytics
SUBSET_TO_ABCAM <- yaml$subset_only_abcam_catalog
STAR_FOR_ABCAM <- yaml$star_for_abcam_product
ONLY_PROTEIN_CODING <- yaml$only_protein_coding

#' Define main entrypoint of app
#' 
#' @export
#' 
#' @import shiny shinytest2
#' @importFrom shinyalert useShinyalert shinyalert
#' @importFrom DT DTOutput renderDT
#' @importFrom tidyr drop_na
#' @importFrom waiter useWaitress Waitress
#' @importFrom cicerone use_cicerone Cicerone
#' @importFrom sortable add_rank_list bucket_list
#' @import fontawesome
#' @import emojifont
#' @import reactable
#' @import htmltools
#' @importFrom rlang is_empty
#' @importFrom shinyanimate withAnim startAnim
#' @importFrom DT datatable
#' @importFrom readr write_lines read_tsv read_csv
#' @importFrom dplyr desc mutate_if distinct
#' @importFrom bsplus use_bs_popover shinyInput_label_embed shiny_iconlink bs_embed_tooltip use_bs_tooltip
#' @importFrom shinyjs useShinyjs hidden toggle reset delay
#' @importFrom grDevices dev.off pdf
#' @importFrom zip zip
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom stringr str_split_fixed
#' @importFrom magrittr set_names
#' @importFrom shinydashboard box dashboardBody dashboardHeader dashboardSidebar
#' dashboardPage menuItem sidebarMenu sidebarMenuOutput tabItem tabItems
#' valueBoxOutput renderMenu updateTabItems tabBox renderValueBox valueBox
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom yaml read_yaml
#' @importFrom feather read_feather
#' @importFrom rdrop2 drop_download
#' @importFrom lubridate with_tz
#' @export cytomarker cytomarker
#' 
#' @param ... Additional arguments
cytomarker <- function(...) {
  
  rm(list = ls())
  
  antibody_info <- read_feather(system.file("catalog.feather", package = "cytomarker")) |>
    mutate(`Protein Expression`  = ifelse(!is.na(ensgene), paste("https://www.proteinatlas.org/", 
                                                                 ensgene,
                                                                 "-", Symbol, "/tissue", sep = ""), NA))
  
  # options(MulticoreParam=quote(MulticoreParam(workers=availableCores())))
  
  compartments <- yaml::read_yaml(system.file("ts_compartments.yml", 
                                              package = "cytomarker"))
  
  
  # applications_parsed <- cytomarker_data$applications
  
  grch38 <- cytomarker_data$grch38
  
  full_palette <- create_global_colour_palette()
  
  markdown_report_path <- system.file(file.path("report", "rmarkdown-report.Rmd"), 
                                      package = "cytomarker")
  
  # devtools will find the file in the inst directory (move to top level)
  cytomarker_token <- readRDS(system.file(file.path("token.rds"),
                                          package = "cytomarker"))
  
  
  dataset_labels <- curated_datasets$tissue |> set_names(gsub("_", " ", 
                                                              curated_datasets$tissue))
  
  other_dataset_labels <- other_curated$tissue |> set_names(gsub("_", " ", 
                                                                 other_curated$tissue))
  
  all_curated_dataset_labels <- c(dataset_labels, other_dataset_labels)
  
  all_zones <- cytomarker_data$time_zones
  # google analytics event tracking: 
  # https://www.gravitatedesign.com/blog/event-tracking-google-analytics/
  # https://absentdata.com/track-external-links/
  ui <- tagList(
    # useWaiter(),
    useWaitress(),
    # attendantBar("progress-bar"),
    includeCSS(system.file(file.path("www", "cytomarker.css"),
                           package = "cytomarker")),
    # only permit the google analytics non-interactively
    if (isTruthy(USE_ANALYTICS) & (!interactive())) {
      tags$head(HTML("<script async src='https://www.googletagmanager.com/gtag/js?id=G-B26X9YQQGT'></script>
          <script>
window.dataLayer = window.dataLayer || [];
function gtag(){dataLayer.push(arguments);}
gtag('js', new Date());
gtag('config', 'G-B26X9YQQGT');

</script>"),
                # tags$script(HTML("window.onbeforeunload = function() {return 'Please visit https://www.surveymonkey.com/ before you leave!';};"))
      )
    },
    # tags$script(HTML("window.onbeforeunload = function() {return 'Please visit https://www.surveymonkey.com/ before you leave!';};"))
    tags$head(tags$style(".modal-dialog{ width:750px}")),
    tags$head(
      tags$script(HTML("window.onbeforeunload = function() {return 'Your changes will be lost!';};"))
    ),
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
    use_cicerone(),
    withAnim(),
    
    dashboardPage(title="cytomarker",
                  
                  # Title
                  dashboardHeader(
                    title = imageOutput("cytomarker_logo"),
                    # span(style = "padding-bottom: 50px;", ), 
                    
                    tags$li(class = "dropdown",
                            hidden(div(id = "session_info",
                                       htmlOutput("current_session_info",
                                                  style = "width:100%; color:white; height:70%; margin-right:6px; margin-top:6.5px;
                                  margin-bottom:-2px; color: white;")))),
                    hidden(tags$li(id = "gene_alias_tag", class = "dropdown",
                                   actionLink("gene_alias_table_viewer", "View gene aliases",
                                              width = "90%"), style = "margin-top:-0.5px;
                                      margin-right: 3px; margin-bottom:-2px; color: white;")),
                    tags$li(class = "dropdown", 
                            actionLink("help_guide", "Guided tour",
                                       width = "90%", style = "margin-top:-0.5px;
                                      margin-left: -3px; margin-right: -10px; margin-bottom:-2px; color: white;")),
                    tags$li(class = "dropdown", 
                            actionLink("time_zone", "Set Time Zone",
                                       width = "90%", style = "margin-top:-0.5px;
                                      margin-right: 15px; margin-bottom:-2px; color: white;"))),
                  
                  dashboardSidebar(
                    sidebarMenu(id = "tabs",
                                # width = 300,
                                menuItem("Get Started", tabName = "inputs", icon = icon("gear")),
                                sidebarMenuOutput(outputId = 'output_menu'),
                                # menuItem("Documentation", tabName = "documentation", icon = icon("bookmark")),
                                column(1, column(11, align = "left", style = "margin-top: 10px", offset = 0, icon("bookmark", style = "margin-left: -10px;"),
                                                 tags$a(href="http://camlab-bioml.github.io/cytomarker-doc/docs/intro",
                                                        " Documentation",
                                                        id = "cytomarker Documentation",
                                                        target="_blank", style = "margin-left: 2px;"))),
                                br(),
                                column(12, align = "center", offset = 0,
                                       hidden(div(id = "analysis_button", actionButton("start_analysis", "Run analysis!",
                                                                                       icon=icon("play", style = "color:black;"),
                                                                                       style = "margin-left: -5px; margin-top: 20px;")))),
                                column(12, align = "center", offset = 0,
                                       hidden(div(id = "download_button", downloadButton("downloadData",
                                                                                         "Save panel", style = "color:black;
                         margin-top: 15px; margin-left:-10px; width: 60%; margin-bottom: 15px;")))),
                                br(),
                                column(1, column(11, align = "left", style = "margin-top: 10px", offset = 0, icon("comments", style = "margin-left: -12px;"),
                                                 tags$a(href="https://www.surveymonkey.com/r/M9X6KPT",
                                                        " Leave Feedback",
                                                        target="_blank", style = "margin-left: 2px;")))
                                
                    )),
                  dashboardBody(
                    # https://stackoverflow.com/questions/52198452/how-to-change-the-background-color-of-the-shiny-dashboard-body
                    tags$head(tags$style(HTML('
                              /* body */
                              .content-wrapper, .right-side {
                              background-color: #FFFFFF;
                              }
                              
                              .skin-blue .main-header .logo {
                              background-color: #FFFFFF; padding-bottom: 10px;
                              }
                              
                              .skin-blue .main-header .logo:hover {
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
                              
                              .small-box {height: 125px}
                              
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
                              fluidRow(column(4, box(id = "input_box",
                                                     title = "Upload scRNA-seq data", status = "primary",
                                                     width = 12, style = "padding-bottom: -15px",
                                                     div(style = "margin-top: -30px; margin-bottom: -25px"),
                                                     fileInput("input_scrnaseq",
                                                               label = "", accept = c(".rds")),
                                                     div(style = "margin-top: -30px")
                              )),
                              tags$h4("or", style = "width = 80%; padding-left: -10px; padding-top: -40px"),
                              div(style = "margin-top: -56px"),
                              column(3, align = "center", style = "margin-left: 15px", box(title = "Select a curated dataset", status = "primary",
                                                                                           width = 12, actionButton("curated_dataset", "Browse human datasets",
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
                              fluidRow(column(4, box(id = "metadata_box", title = "Browse dataset metadata",
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
                                                              numericInput("subset_number", "Number of cells to subsample:", 2000, 
                                                                           min = 100, max = 10000, step = NA, width = NULL) %>%
                                                                shinyInput_label_embed(
                                                                  icon("circle-info") %>%
                                                                    bs_embed_tooltip(title = get_tooltip('subset_number'),
                                                                                     placement = "right")),
                                                              hidden(div(id = "precomputed",
                                                                         checkboxInput("precomputed_dim", "Use precomputed UMAP", value = F))),
                                                              hidden(div(id = "precomputed_list",
                                                                         selectInput("possible_precomputed_dims",
                                                                                     NULL,
                                                                                     choices = NULL,
                                                                                     selected = NULL,
                                                                                     multiple = F))),
                                                              # hidden(div(id = "apps",
                                                              # TODO: for now, applications are blocked because of the catalog format
                                                              hidden(div(selectInput("select_aa", "Antibody applications:", 
                                                                                     NULL, multiple=TRUE) %>%
                                                                           shinyInput_label_embed(
                                                                             icon("circle-info") %>%
                                                                               bs_embed_tooltip(title = get_tooltip('applications'),
                                                                                                placement = "right", html = "true"))))
                                                              # ))
                                                            ),
                                                            bsCollapsePanel(title = HTML(paste0(
                                                              "Previous/custom analysis", tags$span(icon("sort-down",
                                                                                                         style = "position:right; margin-left: 4px; margin-top: -4px;")))), style = "info",
                                                              fileInput("read_back_analysis",
                                                                        label = "Select a yml file from a previous run or upload a custom .txt file of markers",
                                                                        accept = c(".yml", ".txt")) %>%
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
                              fluidRow(column(4, actionButton("markers_change_modal", "Add markers to panel")),
                                       column(8, div(style = "margin-top: -8.5px;"),
                                              radioButtons("panel_sorter", label = "Select marker order", inline = TRUE,
                                                           choices = c("Group by cell type", "Sort alphabetically"),
                                                           selected="Group by cell type"))
                              ),
                              hr(),
                              fluidRow(column(7, hidden(div(id = "marker_visualization",
                                                            tags$span(icon("circle-info")
                                                                      %>%
                                                                        bs_embed_tooltip(title = get_tooltip('marker_visualization'),
                                                                                         placement = "right", html = TRUE)))),),
                                       column(5, hidden(div(id = "marker_display",
                                                            tags$span(icon("circle-info")
                                                                      %>%
                                                                        bs_embed_tooltip(title =  if(isTruthy(STAR_FOR_ABCAM)) get_tooltip('marker_display_abcam') else get_tooltip('marker_display_all'),
                                                                                         placement = "right", html = TRUE)))))),
                              fluidRow(column(4, align = "center", div(style="display: inline-block; font-size: 15px", 
                                                                       htmlOutput("scratch_marker_counts"))),
                                       column(4, align = "center", div(style="display: inline-block; font-size: 15px", 
                                                                       htmlOutput("selected_marker_counts")))),
                              fluidRow(column(8, uiOutput("BL"),
                                              style = "margin-bottom:0px;"
                              ),
                              column(4, align = "left", box(div(style = "margin-left: -40px; overflow-y: scroll;"),
                                                            width = 12,
                                                            height = "1300px",
                                                            status = "primary",
                                                            plotOutput("legend", width = "auto",
                                                                       height = "1300px",
                                                                       
                                                            )),
                                     style = "margin-top:-70px; margin-left: -50px;"))),
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
                                                                                                                      "Gene Marker"),
                                                                                selected = "Cell Type", multiple=FALSE),
                                          column(6, hidden(div(id = "umap_panel_cols",
                                                               selectizeInput("umap_panel_options", 
                                                                              "Color Heatmap by Gene Marker", 
                                                                              choices = NULL, multiple = F,
                                                                              width = "89%",
                                                               )))),
                                          column(2, div(style = "margin-top: 30px; margin-left: 60px;"),
                                                 checkboxInput("show_umap_legend", "Show UMAP plot legends", T,
                                                 ))),
                              fluidRow(column(6, style = "margin-right: -15px; margin-top: 10px",
                                              plotly::plotlyOutput("all_plot", width="550px", height="500px",
                                              )),
                                       column(6, style = "margin-left: -15px; margin-top: 10px",
                                              plotly::plotlyOutput("top_plot", width="550px", height="500px",
                                              )))
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
                              fluidRow(column(12, plotly::plotlyOutput("heatmap", height="auto", width = "auto"))),
                              br()
                      ),
                      
                      tabItem("Metrics",
                              # icon = icon("ruler"),
                              br(),
                              fluidRow(
                                column(10, helpText(get_tooltip('metrics')) %>% tags$span(icon("circle-info") %>%
                                                                                            bs_embed_tooltip(title = get_tooltip('metrics_explanation'),
                                                                                                             placement = "right", html = TRUE)), div(style = "height: 100px;")),
                                column(2,
                                       valueBoxOutput("current_metric_score", width = NULL
                                       ))),
                              # textOutput("cells_per_category"),
                              fluidRow(
                                column(8, plotly::plotlyOutput("metric_plot", height="550px")),
                                column(4, align = "center",
                                       tabBox(width = NULL,
                                              id = "metrics_toggle",
                                              tabPanel(id = "current_metrics_viewed", title = "Current Run Metrics",
                                                       div(style = "align: center;"),
                                                       DTOutput("current_run_metrics", width = "100%")
                                              ),
                                              tabPanel(id = "previous_metrics_viewed", title = "Previous Run Metrics",
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
                              fluidRow(column(12, DTOutput("alternative_markers")))
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
                      tabItem(id = "documentation", "documentation",
                              htmlOutput("cytomarker_hyperlink")
                              # h2("Documentation Preview:"),
                              # htmlOutput("cytomarker_preview"))
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
    heatmap_for_report <- reactiveVal()
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
    
    multimarkers <- reactiveVal(NULL)
    marker_sort <- reactiveVal(NULL)
    
    display <- reactiveVal() # Display options for heatmap
    
    replacements <- reactiveVal() # Where to store table of alternative markers
    
    suggestions <- reactiveVal() # Where to store suggested markers to remove
    
    cells_per_type <- reactiveVal() ## Number of cells per condition/type
    
    summary_cat_tally <- reactiveVal()
    
    cytomarker_palette <- reactiveVal()
    
    proceed_with_analysis <- reactiveVal(TRUE)
    
    specific_cell_types_selected <- reactiveVal()
    
    cell_types_high_enough <- reactiveVal()
    
    cell_types_to_keep <- reactiveVal()
    
    ts_compartments <- reactiveVal(tolower(names(compartments)))
    compartments_selected <- reactiveVal()
    
    ## Current data frame of selected cell type markers and the reactable
    current_cell_type_marker_fm <- NULL
    selected_cell_type_markers <- reactive(getReactableState("cell_type_marker_reactable", "selected"))
    
    any_cells_present <- reactiveVal(TRUE)
    
    num_markers_in_selected <- reactiveVal()
    num_markers_in_scratch <- reactiveVal()
    marker_filtration <- reactiveVal()
    
    cell_types_excluded <- reactiveVal()
    
    marker_suggestions <- reactiveVal()
    alternative_marks <- reactiveVal()
    original_panel <- reactiveVal(NULL)
    
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
    time_zone_set <- reactiveVal(NULL)
    
    # run log variables #
    current_run_log <- reactiveVal()
    previous_run_log <- reactiveVal(NULL)
    previous_run_log_2 <- reactiveVal(NULL)
    current_overall_score <- reactiveVal(NULL)
    previous_overall_score <- reactiveVal(NULL)
    metrics_being_viewed <- reactiveVal(FALSE)
    
    proper_organism <- reactiveVal(TRUE)
    
    downloaded_content <- reactiveVal(FALSE)
    plots_for_markdown <- reactiveVal()
    markers_with_type <- reactiveVal()
    
    curated_selection <- reactiveVal(NULL)
    aliases_table <- reactiveVal(NULL)
    gene_aliases_to_show <- reactiveVal(NULL)
    
    addResourcePath('report', system.file('report', package='cytomarker'))
    
    output$cytomarker_logo <- shiny::renderImage({
      list(src=system.file(file.path("report", "cytomarker-logo.png"), package = "cytomarker"),
           width = "90%",
           height = "9.75%",
           class = "topimg",
           alt = "cytomarker")
    }, deleteFile = F)
    
    output$cytomarker_hyperlink <-  renderUI({
      # url <- a("cytomarker Documentation", href="http://camlab-bioml.github.io/cytomarker-doc/docs/intro")
      # tagList("cytomarker Documentation", url)
      tags$a(href="http://camlab-bioml.github.io/cytomarker-doc/docs/intro",
             "Click Here to access the cytomarker documentation in a new tab",
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
    
    observeEvent(input$help_guide, {
      
      req(input$tabs == "inputs")
      
      guide <- Cicerone$
        new(allow_close = T)$
        step(
          "input_box",
          "Step #1, Option A: Select your own data to analyze",
          "Select a dataset from your local console to load into cytomarker. 
        Use the black information hovertips to get more
      information about what types of data files are compatible.",
        )$
        step(
          "curated_dataset",
          "Step #1, Option B: Select a pre-annotated human dataset",
          "Browse a selection of 24 human tissue types from the Tabula Sapiens Initiative 
        that are pre-annotated and ready for analysis."
        )$
        step(
          "metadata_box",
          "Step #2: Choose a relevant metadata variable to analyze",
          "Browse the metadata of your dataset to find a suitable biological category to analyze, such as
        cell type, disease state, patient status, etc. Use the black information hovertips to 
        see examples of suitable metadata types."
        )$
        step(
          "advanced_collapse",
          "Step #3 (Optional): Customize the run parameters",
          "Modify thresholds for data subsetting, filtering, panel size, and 
        additional run configurations. Additionally, you may upload a previous analysis 
        or a custom list of markers to generate a panel."
        )$
        # step(
        #   "cytomarker_hyperlink",
        #   "Refer to the official documentation",
        #   "Refer to the official cytomarker documentation for additional information on data formatting,
        #   processing, and interpreting outputs."
        # )$
        step(
          "start_analysis",
          "Step #4: Run analysis",
          "Hit <b> Run analysis! </b> to generate an initial panel and 
        populate the dashboard with additional analysis & visualization tools."
        )
      
      
      guide$init(session = session)
      guide$start()
    })
    
    
    observeEvent(input$gene_alias_table_viewer, {
      
      req(allowed_genes())
      req(sce())
      req(aliases_table())
      
      output$table_of_gene_aliases <- DT::renderDataTable(
        gene_aliases_to_show()
      )
      
      showModal(show_alias_table_modal())
      
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
    
    observe(toggle(id = "session_info", condition = isTruthy(time_zone_set())))
    observe(toggle(id = "gene_alias_tag", condition = isTruthy(gene_aliases_to_show())))
    
    # output$cytomarker_preview <- renderUI({
    #   tags$iframe(src="http://camlab-bioml.github.io/cytomarker-doc/docs/intro", height=600, width=1100)
    # })
    
    observeEvent(input$curated_dataset, {
      
      if (!isTruthy(compartments_selected())) compartments_selected(names(compartments)) 
      
      showModal(curated_dataset_modal(names(compartments),
                                      compartments_selected(), names(dataset_labels)))
      
      updateRadioButtons(session, "curated_selection_choice",
                         selected = "Tabula Sapiens")
      
      toggle(id = "curated_ts_compartments", condition = (input$curated_selection_choice == "Tabula Sapiens"))
    })
    
    observeEvent(input$curated_selection_choice, {
      
      toggle(id = "curated_ts_compartments", condition = 
               (input$curated_selection_choice == "Tabula Sapiens"))
      
      if (input$curated_selection_choice == "Other Human") {
        updateSelectInput(session, "curated_options", choices = names(other_dataset_labels),
                          label = "Select a human dataset")
      } else {
        updateSelectInput(session, "curated_options", choices = names(dataset_labels),
                          label = "1. Choose a tissue type to analyze")
      }
      
    })
    
    
    observeEvent(input$curated_compartments, {
      
      possible_curated <- c()
      
      for (sub in input$curated_compartments) {
        possible_curated <- c(possible_curated, compartments[[sub]])
      }
      
      if (length(input$curated_compartments) < 1) possible_curated = names(dataset_labels)      
      updateSelectInput(session, "curated_options", 
                        choices = sort(unique(possible_curated)),
                        selected = ifelse(input$curated_options %in% possible_curated,
                                          input$curated_options, sort(unique(possible_curated))[1]))
      
      ts_compartments(tolower(input$curated_compartments))
      compartments_selected(input$curated_compartments)
      
    }, ignoreNULL = F)
    
    observeEvent(input$curated_options, {
      
      if (input$curated_selection_choice == "Tabula Sapiens") {
        tissue_lab <- dataset_labels[names(dataset_labels) == input$curated_options]
        
        dataset_preview <- paste("<b>", "<a href='https://tabula-sapiens-portal.ds.czbiohub.org/' target='_blank'",
                                 ">", "Tabula Sapiens dataset: ", input$curated_options, "</a>",
                                 "</b>", "<br/>", "<b>", "Cells: ", "</b>", subset(curated_datasets, tissue == tissue_lab)$num_cells[1], 
                                 "<br/>", "<b>", "Genes: ", "</b>", subset(curated_datasets, tissue == tissue_lab)$num_genes[1], "<br/>", 
                                 "<b>", "Cell category of interest: ", "</b>", "cell_ontology_glass",
                                 "<br/>", "<b>", "Cell Type Distribution: ", "</b>", "<br/>", subset(curated_datasets, tissue == tissue_lab)$preview[1],
                                 sep = "")
      } else {
        tissue_lab <- other_dataset_labels[names(other_dataset_labels) == input$curated_options]
        
        dataset_preview <- paste("<b>", "Human dataset: ", input$curated_options, "</a>",
                                 "</b>", "<br/>", "<b>", "Cells: ", "</b>", subset(other_curated, tissue == tissue_lab)$num_cells[1], 
                                 "<br/>", "<b>", "Genes: ", "</b>", subset(other_curated, tissue == tissue_lab)$num_genes[1], "<br/>", 
                                 "<b>", "Cell category of interest: ", "</b>", subset(other_curated, tissue == tissue_lab)$category_interest[1],
                                 "<br/>", "<b>", "Cell Type Distribution: ", "</b>", "<br/>", subset(other_curated, tissue == tissue_lab)$preview[1],
                                 sep = "")
      }
      
      output$curated_set_preview <- renderPrint({HTML(dataset_preview)})
      
    })
    
    observeEvent(input$pick_curated, {
      req(input$curated_options)
      req(input$curated_selection_choice)
      
      # future_promise({
      removeModal()
      reset("input_scrnaseq")
      pre_upload_configuration()
      umap_all(NULL)
      umap_top(NULL)
      
      tissue_lab <- all_curated_dataset_labels[names(
        all_curated_dataset_labels) == input$curated_options]
      
      curated_selection(input$curated_options)
      
      withProgress(message = 'Configuring curated selection', value = 0, {
        setProgress(value = 0)
        
        if (!file.exists(file.path(tempdir(), paste(tissue_lab, ".rds", sep = "")))) {
          incProgress(detail = "Downloading curated dataset")
          origin_dir <- if (input$curated_selection_choice == "Tabula Sapiens")
            "tabula_sapiens/" else "other_curated_cytomarker/"
          
          drop_download(paste(origin_dir, 
                              paste(tissue_lab, ".rds", sep = ""),
                              sep = ""),
                        local_path = tempdir(),
                        overwrite = T,
                        dtoken = cytomarker_token)
        }
        
        setProgress(value = 0.5)
        incProgress(detail = "Reading input dataset")
        
        input_sce <- read_input_scrnaseq(file.path(tempdir(), 
                                                   paste(tissue_lab, ".rds", sep = "")))
        
        
        if (input$curated_selection_choice == "Tabula Sapiens") {
          if (length(ts_compartments()) < 1) ts_compartments(tolower(names(compartments)))
          
          input_sce <- input_sce[,input_sce$compartment %in% ts_compartments()]
          input_sce$compartment <- factor(input_sce$compartment,
                                          levels = ts_compartments())
          default_category_curated("cell_ontology_class")
        } else {
          default_category_curated(subset(other_curated, tissue == tissue_lab)$category_interest[1])
        }
        incProgress(detail = "Parsing gene names and assays")
        
        setProgress(value = 1)
      })
      
      post_upload_configuration(input_sce)
      
      update_metadata_column()
      
      startAnim(session, 'analysis_button', 'rubberBand')
      
    })
    # })
    
    
    ### UPLOAD FILE ###
    observeEvent(input$input_scrnaseq, {
      
      # future_promise({
      default_category_curated(NULL)
      curated_selection(NULL)
      umap_all(NULL)
      umap_top(NULL)
      withProgress(message = 'Configuring input selection', value = 0, {
        setProgress(value = 0)
        pre_upload_configuration()
        setProgress(value = 0.25)
        incProgress(detail = "Reading input dataset")
        input_sce <- read_input_scrnaseq(input$input_scrnaseq$datapath)
        if (!isTruthy(input_sce)) {
          invalid_modal()
          proceed_with_analysis(FALSE)
        }
        
        req(proceed_with_analysis())
        
        req(isTruthy(input_sce))
        incProgress(detail = "Parsing gene names and assays")
        setProgress(value = 0.85)
        found_human <- check_for_human_genes(input_sce)
        if (!isTruthy(found_human)) {
          throw_error_or_warning(type = "error", message = "cytomarker has detected that the gene names 
                             provided are not human. Currently, only human datasets are supported.")
          proper_organism(found_human)
        }
        setProgress(value = 1)
      })
      
      req(proper_organism())
      post_upload_configuration(input_sce)
      
      startAnim(session, 'analysis_button', 'rubberBand')
      
    })
    # })
    
    to_listen_preview <- reactive({
      list(input$pick_curated, input$coldata_column, input$input_scrnaseq)
    })
    
    observeEvent(to_listen_preview(), {
      req(input$coldata_column)
      
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
    })
    
    
    observeEvent(input$coldata_column, {
      req(sce())
      
      update_metadata_column()
      
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
    
    observeEvent(input$select_aa, {
      
      # assume that selecting the applications creates a valid panel for analysis
      valid_existing_panel(TRUE)
      set_allowed_genes()
      if (isTruthy(aliases_table())) {
        gene_aliases_to_show(aliases_table()[aliases_table()$`Alias in dataset` %in% allowed_genes(),])
      }
      
    }, ignoreNULL = F)
    
    observeEvent(input$show_cat_table, {
      req(sce())
      req(pref_assay())
      req(input$coldata_column)
      
      cell_category_table <- create_table_of_hetero_cat(sce(),
                                                        input$coldata_column)
      
      # if (!isTruthy(cell_min_threshold())) cell_min_threshold(2)
      
      
      # if (!isTruthy(specific_cell_types_selected())) specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      
      
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
      
      showModal(markers_add_modal(allowed_genes(), names(fms()[[1]]), session))
    })
    
    # ### ANTIBODY EXPLORER ###
    output$antibody_table <- renderReactable({
      req(current_markers())
      req(df_antibody())
      
      
      reactable(df_antibody(),
                searchable = TRUE,
                filterable = TRUE,
                groupBy = "Symbol",
                columns = list(
                  `Vendor` = colDef(aggregate = "unique",
                                                       filterInput = function(values, name) {
                                                         tags$select(
                                                           # Set to undefined to clear the filter
                                                           onchange = sprintf("Reactable.setFilter('antibody-select', '%s', event.target.value || undefined)", name),
                                                           # "All" has an empty value to clear the filter, and is the default option
                                                           tags$option(value = "", "All"),
                                                           tags$html(multiple = T),
                                                           lapply(unique(values), tags$option),
                                                           "aria-label" = sprintf("Filter %s", name),
                                                           style = "width: 100%; height: 28px"
                                                         )
                                                       }
                ),
                # `Target` = colDef(aggregate = "unique"),
                #                 `Listed Applications` = colDef(
                #                   filterable = TRUE,
                #                   # Filter by case-sensitive text match
                #                   filterMethod = JS("export function MultiSelectFilterFn(rows, id, filterValues) {
                #     if (filterValues.length === 0) return rows;
                # 
                #     return rows.filter(r => filterValues.includes(r.values[id]));
                # }")
                #                 ),
                `Human Protein Atlas` = colDef(html = T)
                ),
                sortable = TRUE,
                elementId = "antibody-select")
    })
    
    ### PLOTS ###
    output$all_plot <- plotly::renderPlotly({
      # req(umap_all())
      req(column())
      req(plots$all_plot)
      
      columns <- column()
      
      plots$all_plot
    })
    
    output$top_plot <- plotly::renderPlotly({
      # req(umap_top())
      req(column())
      req(plots$top_plot)
      
      plots$top_plot
      
    })
    
    output$heatmap <- plotly::renderPlotly({
      req(heatmap())
      req(column())
      
      heatmap()
    })
    
    output$metric_plot <- plotly::renderPlotly({
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
    
    #### Overall Score tables ######
    
    to_listen_metrics <- reactive({
      list(input$tabs, input$metrics_toggle)
    })
    
    
    observeEvent(to_listen_metrics(),{
      req(current_overall_score())
      
      metrics_val <- ifelse(input$metrics_toggle == "Current Run Metrics", current_overall_score()$score,
                            ifelse(!is.null(previous_overall_score()), previous_overall_score()$score,
                                   "None"))
      
      metrics_counts <- ifelse(input$metrics_toggle == "Current Run Metrics", current_overall_score()$counts,
                               ifelse(!is.null(previous_overall_score()), previous_overall_score()$counts,
                                      "None"))
      
      which_set <- ifelse(input$metrics_toggle == "Current Run Metrics", "Current Run, ", "Previous Run, ")
      
      if (input$tabs == "Metrics") {
        
        output$current_metric_score <- renderValueBox({
          valueBox(
            subtitle = paste0("Accuracy, ", which_set, "n = ", metrics_counts),
            value = metrics_val,
            icon = NULL,
            width = NULL,
            color = "blue",
          )
        })
        
      }
      
    })
    
    ### ANALYSIS ###
    observeEvent(input$start_analysis, {
      req(input$panel_size)
      req(input$coldata_column)
      req(sce())
      
      showNotification("Starting analysis",
                       type = 'message',
                       duration = 2)
      
      library(caret, quiet = T)
      library(Seurat, quiet = T)
      library(clustifyr, quiet = T)
      library(ggplot2, quiet = T)
      library(plotly, quiet = T)
      
      proceed_with_analysis(TRUE)
      
      if (isTruthy(input$precomputed_dim)) {
        use_precomputed_umap(TRUE)
        umap_precomputed_col(input$possible_precomputed_dims)
      }
      
      if (input$panel_size < 2) {
        panel_too_small_modal(input$panel_size)
        proceed_with_analysis(FALSE)
      }
      
      req(proceed_with_analysis())
      
      if (!is_empty(input$bl_top)) {
        current_panel_not_valid <- setdiff(input$bl_top, rownames(sce()))
        if (length(current_panel_not_valid) > 0) {
          valid_existing_panel(FALSE)
          current_pan_not_valid_modal(current_panel_not_valid, "current dataset:")
        }
        
        use_precomputed_umap(input$precomputed_dim)
        
        set_allowed_genes()
        
        current_panel_not_allowed <- setdiff(input$bl_top, allowed_genes())
        
        if (length(current_panel_not_allowed) > 0) {
          valid_existing_panel(FALSE)
          current_pan_not_valid_modal(current_panel_not_allowed, "gene lists for the chosen antibody applications:")
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
      
      if (nrow(cell_cat_summary) > 100) {
        proceed_with_analysis(FALSE)
        invalid_metadata_modal(input$coldata_column)
      } else {
        proceed_with_analysis(TRUE)
      }
      
      req(proceed_with_analysis())
      
      cell_types_high_enough(remove_cell_types_by_min_counts(cell_cat_summary,
                                                             sce(), 
                                                             input$coldata_column, 
                                                             cell_min_threshold()))
      
      cell_types_to_keep(intersect(specific_cell_types_selected(),
                                   cell_types_high_enough()))
      
      sce(remove_null_and_va_from_cell_cat(sce(), input$coldata_column))
      
      # showNotification("Ignoring any cells with input category set to NA or null",
      #                  type = 'message',
      #                  duration = 4)
      
      sce(create_sce_column_for_analysis(sce(), cell_types_to_keep(), 
                                         input$coldata_column))
      
      if (dim(sce()[,sce()$keep_for_analysis == "Yes"])[2] == 0) {
        no_cells_left_modal()
        any_cells_present(FALSE)
      } else {
        any_cells_present(TRUE)
      }
      
      req(proceed_with_analysis())
      req(any_cells_present())
      
      ## Set initial markers:
      scratch_markers_to_keep <- input$bl_scratch
      
      columns <- good_col(sce()[,sce()$keep_for_analysis == "Yes"], input$coldata_column)
      column(columns$good)
      col <- columns$bad
      
      
      if (input$subsample_sce) {
        # incProgress(detail = "Subsampling data")
        
        avail_cols <- which(sce()$keep_for_analysis == "Yes")
        to_subsample <- sample(avail_cols, min(length(avail_cols), min(ncol(sce()), input$subset_number)),
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
      
      req(proceed_with_analysis())
      
      new_waitress <- Waitress$new(theme = "overlay-percent", infinite = TRUE)
      new_waitress$start()
      
      different_cells <- setdiff(specific_cell_types_selected(),
                                 cell_types_high_enough())
      
      if (length(different_cells) > 0 & isTruthy(any_cells_present())) {
        output$ignored_cell_types <- renderDataTable(data.frame(
          `Cell Type Ignored` = different_cells))
        cell_types_excluded(different_cells[!is.null(different_cells)])
        cell_type_ignored_modal(cell_min_threshold(),
                                cell_types_excluded())
      }
      
      
      # setProgress(value = 1)
      
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
          
          set_allowed_genes()
          
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
            markers <- list(recommended_markers = input$bl_recommended[input$bl_recommended %in% allowed_genes()],
                            scratch_markers = input$bl_scratch[input$bl_scratch %in% allowed_genes()],
                            top_markers = input$bl_top[input$bl_top %in% allowed_genes()])
            
          } else {
            ## We compute the set of markers for the first time
            
            if (isTruthy(markers_reupload())) {
              
              markers <- list(recommended_markers = 
                                markers_reupload()$top_markers[markers_reupload()$top_markers %in% allowed_genes()],
                              top_markers = markers_reupload()$top_markers[markers_reupload()$top_markers 
                                                                           %in% allowed_genes()],
                              scratch_markers = markers_reupload()$scratch_markers[markers_reupload()$scratch_markers 
                                                                                   %in% allowed_genes()])
              
              
            } else {
              markers_list <- get_markers(fms(), 
                                          # Adding 10 to make sure panel size is approximate 
                                          # since a) the same marker is selected multiple times and
                                          # b) excess markers are removed
                                          input$panel_size, 
                                          input$marker_strategy, 
                                          sce()[,sce()$keep_for_analysis == "Yes"],
                                          allowed_genes())
              
              markers <- markers_list$marker[!is.na(markers_list$marker)]
              markers$scratch_markers <- scratch_markers_to_keep
              multimarkers(markers_list$multimarkers)
              
              if (isTruthy(markers_list$missing)) {
                throw_error_or_warning(type = 'notification',
                                       message = paste("No markers were found for the following cell types: ",
                                                       paste(markers_list$missing, 
                                                             collapse = ", "), 
                                                       ". This is likely because there are too few cells of these types."),
                                       duration = 10)
                cell_types_missing_markers(markers_list$missing)
                
              }
            }
          }
          
          setProgress(value = 1)
          
          if(length(markers$recommended_markers) < input$panel_size &
             (!isTruthy(reupload_analysis()))) {
            showNotification("cytomarker found genes that are good markers for multiple cell types.
            This will result in a smaller panel size than requested.
                               You may manually add additional markers.",
                             type = 'message',
                             duration = NULL)
          }
          current_markers(set_current_markers_safely(markers, fms()))
          
          # SMH
          
          
          if (!isTruthy(original_panel())) {
            original_panel(current_markers()$recommended_markers)
          }
          
          num_markers_in_selected(length(current_markers()$top_markers))
          num_markers_in_scratch(length(current_markers()$scratch_markers))
          cells_per_type(table(SummarizedExperiment::colData(
            sce()[,sce()$keep_for_analysis == "Yes"])[[column()]]))
          
          if (isTruthy(input$select_aa)) {
            clean_table <- antibody_info |> mutate(limit_applications = sapply(antibody_info$`Listed Applications`, 
                                                                               FUN = function(x) in_row(input$select_aa, x))) |>
              dplyr::filter(limit_applications == isTruthy(limit_applications)) |> select(-limit_applications)
          } else {
            clean_table <- antibody_info
          }
          
          df_antibody(dplyr::filter(clean_table, Symbol %in% 
                                      current_markers()$top_markers) |>
                        # dplyr::distinct(ID, .keep_all = T) |>
                        mutate(
                          `Vendor` = factor(`Vendor`),
                               Symbol = factor(Symbol),
                               `Human Protein Atlas` = ifelse(!is.na(`Protein Expression`),
                                                              paste0('<a href="',`Protein Expression`, '"', 'id=', '"', `Product Name`, '"',
                                                                     'onClick=”_gaq.push([‘_trackEvent’, ‘abcam_link’, ‘click’, ‘abcam_link’, ‘abcam_link’]);”',
                                                                     ' target="_blank" rel="noopener noreferrer"',
                                                                     '>', "View on ",
                                                                     as.character(icon("external-link-alt")),
                                                                     "Human Protein Atlas",
                                                                     '</a>'), "None")
                               ) |> dplyr::select(-c(`Protein Expression`, `ensgene`)))
          
          update_analysis()
          updateSelectizeInput(
            session = session,
            inputId = "umap_panel_options",
            # choices = current_markers()$top_markers,
            # selected = current_markers()$top_markers[1])
            choices = allowed_genes(),
            selected = current_markers()$top_markers[1],
            server = T)
          
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
      
      # reset downloaded status when analysis is rerun
      downloaded_content(FALSE)
      toggle(id = "download_button", condition = isTruthy(current_markers()))
      
      umap_all(NULL)
      umap_top(NULL)
      reset_panel(FALSE)
      
      aliases_table(aliases_table()[aliases_table()$`Alias in dataset` %in% 
                                      current_markers()$top_markers | aliases_table()$`Alias in dataset` %in% 
                                      current_markers()$scratch_markers,])
      
      gene_aliases_to_show(gene_aliases_to_show()[gene_aliases_to_show()$Alias %in% 
                                                    current_markers()$top_markers | gene_aliases_to_show()$Alias %in% 
                                                    current_markers()$scratch_markers,])
      
      toggle(id = "gene_alias_tag", condition = isTruthy(gene_aliases_to_show()))
      
      new_waitress$close()
      
      
      
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
         (input$add_markers %in% allowed_genes())) {
        
        # && (! input$add_markers %in% 
        #    current_markers()$top_markers) && (! input$add_markers %in% 
        #                                       current_markers()$scratch_markers)) {
        ## Need to update markers:
        new_marker <- input$add_markers[!input$add_markers %in% current_markers()$top_markers &
                                          !input$add_markers %in% current_markers()$scratch_markers]
        
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
        
      } else if(!(input$add_markers %in% rownames(sce()))) {
        dne_modal(dne = input$add_markers)
      }
    })
    
    
    #### Add genes from expression tab #####
    
    observeEvent(input$add_violin_genes, {
      req(input$genes_for_violin)
      
      if(!is.null(input$genes_for_violin) && 
         all(input$genes_for_violin %in% allowed_genes())) {
        
        to_add <- input$genes_for_violin[!input$genes_for_violin %in% current_markers()$top_markers &
                                           !input$genes_for_violin %in% current_markers()$scratch_markers]
        
        cm <- current_markers()
        markers <- list(recommended_markers = cm$recommended_markers,
                        scratch_markers = input$bl_scratch,
                        top_markers = unique(c(to_add,
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
        updateTabsetPanel(session, "tabs", "marker_selection")
      }
      
    })
    
    observeEvent(input$add_to_selected, { # Add uploaded markers
      req(input$uploadMarkers)
      req(allowed_genes())
      
      uploaded_markers <- readLines(input$uploadMarkers$datapath)
      alias <- get_gene_aliases(uploaded_markers, cytomarker_data$gene_mapping,
                                allowed_genes())
      
      aliases_table(if (isTruthy(aliases_table())) rbind(aliases_table(), alias$table) else alias$table)
      
      gene_aliases_to_show(aliases_table()[aliases_table()$`Alias in dataset` %in% allowed_genes(),])
      
      output$table_of_gene_aliases <- DT::renderDataTable(
        gene_aliases_to_show()
      )
      
      to_add <- alias$merged
      to_add <- to_add[!to_add %in% current_markers()$top_markers &
                         !to_add %in% current_markers()$scratch_markers]
      
      not_sce <- setdiff(uploaded_markers, to_add)
      
      if (isTruthy(uploaded_markers)) {
        markers <- list(recommended_markers = current_markers()$recommended_markers,
                        scratch_markers = input$bl_scratch,
                        top_markers = unique(c(current_markers()$top_markers, to_add)))
        
        current_markers(
          set_current_markers_safely(markers, fms())
        )
        
        num_markers_in_selected(length(current_markers()$top_markers))
        num_markers_in_scratch(length(current_markers()$scratch_markers))
        
        update_BL(current_markers(), num_markers_in_selected(),
                  num_markers_in_scratch(),
                  names(fms()[[1]]))
        
        if(length(not_sce) > 0) {
          # warning_modal(not_sce)
          showNotification(paste("The following genes were not found in the dataset and were not added:",
                                 paste(not_sce, collapse = ",")), duration = 3)
        }
        
      }
      
      if (length(alias$merged) > length(uploaded_markers) | !all(alias$merged %in% uploaded_markers) |
          !all(uploaded_markers %in% alias$merged)) {
        showModal(additional_aliases_modal())
      }
    })
    
    observeEvent(input$replace_selected, { # Replace selected markers by uploaded markers
      req(input$uploadMarkers)
      req(allowed_genes())
      req(current_markers())
      
      uploaded_markers <- readLines(input$uploadMarkers$datapath)
      marker <- uploaded_markers[!uploaded_markers %in% current_markers()$scratch_markers]
      not_sce <- setdiff(uploaded_markers, marker)
      
      alias <- get_gene_aliases(marker, cytomarker_data$gene_mapping,
                                allowed_genes())
      
      aliases_table(alias$table)
      gene_aliases_to_show(aliases_table()[aliases_table()$`Alias in dataset` %in% allowed_genes(),])
      
      output$table_of_gene_aliases <- DT::renderDataTable(
        gene_aliases_to_show()
      )
      
      marker <- alias$merged
      
      if (isTruthy(marker)) {
        markers <- list(recommended_markers = current_markers()$recommended_markers,
                        scratch_markers = input$bl_scratch,
                        top_markers = marker)
        
        current_markers(
          set_current_markers_safely(markers, fms())
        )
        
        if(length(not_sce) > 0) {
          # warning_modal(not_sce)
          showNotification(paste("The following genes were not found in the dataset and were not added:",
                                 paste(not_sce, collapse = ",")), duration = 3)
        }
        
        num_markers_in_selected(length(current_markers()$top_markers))
        num_markers_in_scratch(length(current_markers()$scratch_markers))
        
        update_BL(current_markers(), num_markers_in_selected(),
                  num_markers_in_scratch(),
                  names(fms()[[1]]))
      }
      
      if (length(alias$merged) > length(marker) | !all(alias$merged %in% marker) |
          !all(marker %in% alias$merged)) {
        showModal(additional_aliases_modal())
      }
      
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
      req(fms())
      
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
                                                   !allowed_genes() %in% current_markers()$scratch_markers],
                                 fms()) |>
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
      
      toggle(id = "precomputed_list", condition = isTruthy(input$precomputed_dim))
      updateSelectInput(session, "possible_precomputed_dims", 
                        choices = possible_umap_dims(),
                        selected = possible_umap_dims()[1])
      
      umap_all(NULL)
      
    }, ignoreNULL = F)
    
    
    observeEvent(input$possible_precomputed_dims, {
      req(sce())
      req(input$possible_precomputed_dims)
      
      if (isTruthy(input$precomputed_dim)) {
        use_precomputed_umap(TRUE)
        umap_precomputed_col(input$possible_precomputed_dims)
      }
      
    }, ignoreNULL = F)
    
    ## Re-upload previous analysis
    
    observeEvent(input$read_back_analysis, {
      
      if(!isTruthy(sce())) {
        showNotification("Please upload or select a curated dataset prior to uploading a prior configuration.",
                         duration = 3)
      }
      
      req(sce())
      req(allowed_genes())
      
      if (file_ext(input$read_back_analysis$datapath) == "txt") {
        uploaded_markers <- readLines(input$read_back_analysis$datapath)
        alias <- get_gene_aliases(uploaded_markers, cytomarker_data$gene_mapping,
                                  allowed_genes())
        
        aliases_table(alias$table)
        gene_aliases_to_show(aliases_table()[aliases_table()$`Alias in dataset` %in% allowed_genes(),])
        
        output$table_of_gene_aliases <- DT::renderDataTable(
          gene_aliases_to_show()
        )
        
        if (length(alias$merged) > length(uploaded_markers) | !all(alias$merged %in% uploaded_markers) |
            !all(uploaded_markers %in% alias$merged)) {
          showModal(additional_aliases_modal())
        }
        
        markers_reupload(list(top_markers = alias$merged,
                              scratch_markers = NULL))
      }
      
      
      if (file_ext(input$read_back_analysis$datapath) == "yml") {
        
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
        if (isTruthy(yaml_back$`Target panel size`)) {
          if (yaml_length != yaml_back$`Target panel size` & yaml_length > 0) {
            updateNumericInput(session, "panel_size", value = yaml_length)
          }
        }
        
        if (isTruthy(sce())) {
          if ((isTruthy(yaml_back$`Number of columns (cells)`) &
               isTruthy(yaml_back$`Number of rows (features)`))) {
            if ((yaml_back$`Number of columns (cells)` != ncol(sce()) |
                 yaml_back$`Number of rows (features)` != nrow(sce()))) {
              reupload_warning_modal("Incompatible dimensions", "Number of cells and genes")
            }
          }
          if (isTruthy(yaml_back$`Heterogeneity source`)) {
            if (yaml_back$`Heterogeneity source` %in% colnames(SummarizedExperiment::colData(sce()))) {
              
              
              updateSelectInput(session, "coldata_column", choices = colnames(SummarizedExperiment::colData(sce())),
                                selected = yaml_back$`Heterogeneity source`)
              updateSelectInput(session, "user_selected_cells", 
                                yaml_back$`Cell Types Analyzed`)
              
              types_to_add <- if(all(yaml_back$`Cell Types Analyzed` %in% unique(sce()[[yaml_back$`Heterogeneity source`]]))) 
                yaml_back$`Cell Types Analyzed` else unique(sce()[[yaml_back$`Heterogeneity source`]])
              
              specific_cell_types_selected(types_to_add)
              # cell_types_to_keep(yaml_back$`User selected cells`)
              reupload_cell_types(TRUE)
              
            } else {
              if (isTruthy(yaml_back$`Heterogeneity source`)) {
                reupload_warning_modal("Cell type not found", yaml_back$`Heterogeneity source`)
                specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
                
              }
            }
          } else {
            specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
          }
          
          
          if (isTruthy(yaml_back$`Pre-computed UMAP`)) {
            updateCheckboxInput(session, inputId = "precomputed_dim",
                                value = yaml_back$`Pre-computed UMAP`)
            use_precomputed_umap(yaml_back$`Pre-computed UMAP`)
            
          }
          
          if (isTruthy(use_precomputed_umap())) {
            possible_umap_dims(detect_umap_dims_in_sce(sce()))
            updateSelectInput(session, "possible_precomputed_dims", 
                              choices = possible_umap_dims(),
                              selected = possible_umap_dims()[1])
          }
        }
        if (isTruthy(yaml_back$`Assay used`)) {
          pref_assay(yaml_back$`Assay used`)
        }
        
        if (isTruthy(yaml_back$`Min Cell Category cutoff`)) {
          cell_min_threshold(yaml_back$`Min Cell Category cutoff`)
          updateNumericInput(session, "min_category_count", value = yaml_back$`Min Cell Category cutoff`)
        }
        
        if (isTruthy(yaml_back$`Subsampling Used`)) {
          updateCheckboxInput(session, "subsample_sce", value = yaml_back$`Subsampling Used`)
        }  
        
        if (isTruthy(sce()) & isTruthy(yaml_back$`Subsampling Used`) & 
            isTruthy(yaml_back$`Subsampling number`)) {
          if (yaml_back$`Subsampling number` <= ncol(sce())) {
            updateNumericInput(session, "subset_number", value = yaml_back$`Subsampling number`)
          }
        }
        
        if (isTruthy(yaml_back$`Marker strategy`)) {
          updateRadioButtons(session, "marker_strategy", selected = yaml_back$`Marker strategy`)
        }
        
        if (isTruthy(yaml_back$`Antibody applications`)) {
          updateSelectInput(session, "select_aa", selected = yaml_back$`Antibody applications`)
        }
        
      }
      
      reupload_analysis(TRUE)
      
    })
    
    # observeEvent(input$accept_aliases, {
    #   req(custom_marker_upload())
    #   
    #   showNotification("Updating the current panel with uploaded markers and associated aliases.",
    #                    type = "message", duration = 3)
    #   
    #   markers_reupload(list(top_markers = custom_marker_upload()$merged,
    #                                          scratch_markers = NULL))
    #   removeModal()
    #   
    # })
    # 
    # observeEvent(input$decline_aliases, {
    #   req(custom_marker_upload())
    #   
    #   showNotification("Updating the current panel with uploaded markers",
    #                    type = "message", duration = 3)
    #   
    #   markers_reupload(list(top_markers = custom_marker_upload()$original,
    #                         scratch_markers = NULL))
    #   removeModal()
    #   
    # })
    
    observeEvent(input$create_reset, {
      req(sce())
      req(current_markers())
      
      showModal(reset_analysis_modal())
    })
    
    observeEvent(input$dismiss_marker_reset, {
      proceed_with_analysis(TRUE)
      removeModal()
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
      reset("read_back_analysis")
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                names(fms()[[1]]))
      valid_existing_panel(TRUE)
      original_panel(NULL)
      umap_top(NULL)
      umap_all(NULL)
      removeModal()
      updateTabsetPanel(session, "tabs", "marker_selection")
      proceed_with_analysis(TRUE)
      
    })
    
    observeEvent(input$reset_marker_panel_reupload, {
      showModal(reset_analysis_modal())
      
    })
    
    to_listen_umap <- reactive({
      list(input$umap_options, input$umap_panel_options, input$show_umap_legend,
           input$tabs)
    })
    
    observeEvent(input$tabs, {
      
      req(sce())
      req(current_markers())
      req(column())
      req(!isTruthy(reset_panel()))
      
      if (input$tabs == "UMAP") {
        
        if (!isTruthy(umap_all()) | !isTruthy(umap_top())) {
          umap_waitress <- Waitress$new(theme = "overlay-percent", infinite = TRUE)
          umap_waitress$start()
          
          if (!isTruthy(umap_all())) {
            umap_all(get_umap(sce()[,sce()$keep_for_analysis == "Yes"],
                              column(), pref_assay(), 
                              use_precomputed_umap(),
                              umap_precomputed_col(), F, allowed_genes()))
            
            plots$all_plot <- suppressWarnings(plotly::plot_ly(umap_all(), x=~UMAP_1, y=~UMAP_2, color=~get(input$coldata_column), 
                                                               text=~get(input$coldata_column), 
                                                               type='scatter', hoverinfo="text", colors=cytomarker_palette()) %>% 
                                                 layout(title = "UMAP all genes"))
          }
          
          if (!isTruthy(umap_top())) {
            
            umap_top(get_umap(sce()[,sce()$keep_for_analysis == "Yes"][current_markers()$top_markers,], 
                              column(), pref_assay(), use_precomputed_umap(),
                              umap_precomputed_col(),
                              T, current_markers()$top_markers))
            
            plots$top_plot <- suppressWarnings(plotly::plot_ly(umap_top(), x=~UMAP_1, y=~UMAP_2, color=~get(input$coldata_column), 
                                                               text=~get(input$coldata_column), 
                                                               type='scatter', hoverinfo="text", colors=cytomarker_palette()) %>% 
                                                 layout(title = "UMAP selected markers"))
            
          }
          
          
          
          umap_waitress$close()
          
        }
        
        plots_for_markdown(list(all = suppressWarnings(plotly::plot_ly(umap_all(), x=~UMAP_1, y=~UMAP_2, color=~get(input$coldata_column), 
                                                                       text=~get(input$coldata_column), 
                                                                       type='scatter', hoverinfo="text", colors=cytomarker_palette()) %>% 
                                                         layout(title = "UMAP all genes")),
                                top = plots$top_plot,
                                metric = plots_for_markdown()$metric))
      }
      
    })
    
    
    observeEvent(to_listen_umap(), {
      req(sce())
      req(input$umap_options)
      req(input$umap_panel_options)
      req(umap_all())
      req(umap_top())
      req(column() == input$coldata_column)
      req(current_markers())
      req(!isTruthy(reset_panel()))
      req(input$tabs == "UMAP")
      
      umap_colouring(input$umap_options)
      
      toggle(id = "umap_panel_cols", condition = input$umap_options != "Cell Type")
      
      if (umap_colouring() == "Cell Type") {
        plots$all_plot <- suppressWarnings(plotly::plot_ly(umap_all(), x=~UMAP_1, y=~UMAP_2, color=~get(input$coldata_column),
                                                           text=~get(input$coldata_column), 
                                                           type='scatter', hoverinfo="text", colors=cytomarker_palette()) %>% 
                                             layout(title = "UMAP all genes",
                                                    showlegend = input$show_umap_legend))
        
        plots$top_plot <- suppressWarnings(plotly::plot_ly(umap_top(), x=~UMAP_1, y=~UMAP_2, color=~get(input$coldata_column),
                                                           text=~get(input$coldata_column),
                                                           type='scatter', hoverinfo="text", colors=cytomarker_palette()) %>%
                                             layout(title = "UMAP selected markers",
                                                    showlegend = input$show_umap_legend))
        umap_top_gene(FALSE)
        umap_all_gene(FALSE)
        
      } else {
        
        gene_plot_top <- clustifyr::plot_gene(assay(sce()[,sce()$keep_for_analysis == "Yes"], pref_assay()),
                                              umap_top() |> select(UMAP_1, UMAP_2) |> drop_na(),
                                              input$umap_panel_options)
        
        gene_plot_all <- clustifyr::plot_gene(assay(sce()[,sce()$keep_for_analysis == "Yes"], pref_assay()),
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
        
        plots$top_plot <- suppressWarnings(plotly::plot_ly(umap_top_gene(),
                                                           x = ~UMAP_1, y = ~UMAP_2, 
                                                           color = ~Expression,
                                                           type='scatter',
                                                           text = ~lab,
                                                           hoverinfo = "text",
                                                           colors = c("grey60",
                                                                      "blue")) %>%
                                             # cytomarker_palette()[current_markers()$associated_cell_types[input$umap_panel_options]])) %>%
                                             layout(title = "UMAP selected markers"))
        
        plots$all_plot <- plotly::hide_guides(suppressWarnings(plotly::plot_ly(umap_all_gene(),
                                                                               x = ~UMAP_1, y = ~UMAP_2, 
                                                                               color = ~Expression,
                                                                               type='scatter',
                                                                               text = ~lab,
                                                                               hoverinfo = "text",
                                                                               colors = c("grey60",
                                                                                          "blue")) %>%
                                                                 # cytomarker_palette()[current_markers()$associated_cell_types[input$umap_panel_options]])) %>%
                                                                 layout(title = "UMAP all genes",
                                                                        showlegend = F)))
        
        
      }
      
    }, ignoreNULL = F)
    
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
            suppressWarnings(ggplot2::ggplot(viol_frame, aes(x = `Cell Type`, y = Expression, 
                                                             fill = `Cell Type`)) + ggplot2::geom_violin(alpha = 1) +
                               ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
                               ggplot2::facet_wrap(~Gene, scales="free_y") + 
                               ggplot2::scale_fill_manual(values = cytomarker_palette()))
          })
        } else {
          output$expression_violin <- renderPlot({
            suppressWarnings(ggplot2::ggplot(viol_frame, aes(x = `Cell Type`, y = Expression, 
                                                             fill = Gene)) + ggplot2::geom_violin(alpha = 1) +
                               ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
                               # facet_wrap(~Gene) + 
                               ggplot2::scale_fill_manual(values = full_palette))
          })
          
        }
        
      } else {
        output$expression_violin <- renderPlot({
          ggplot2::ggplot() + ggplot2::theme_void()
        })
      }
    }, ignoreNULL = FALSE)
    
    observeEvent(input$panel_sorter, {
      req(current_markers())
      req(fms())
      
      update_BL(current_markers(), num_markers_in_selected(),
                num_markers_in_scratch(),
                names(fms()[[1]]))
      
      marker_sort(input$panel_sorter)
      
    })
    
    # update the selected metadata column and subtypes
    
    update_metadata_column <- function() {
      column(input$coldata_column)
      
      if (!isTruthy(reupload_cell_types())) {
        updateSelectInput(session, "user_selected_cells",
                          unique(sce()[[input$coldata_column]]))
        specific_cell_types_selected(unique(sce()[[input$coldata_column]]))
      }
      
      cell_types_to_keep(NULL)
      
      reupload_analysis(FALSE)
      reupload_cell_types(FALSE)
    }
    
    
    ### UPDATE SELECTED MARKERS ###
    update_BL <- function(markers, top_size, scratch_size, unique_cell_types,
                          adding_alternative = F) {
      set.seed(12345L)
      
      unique_cell_types <- sort(unique_cell_types)
      
      palette_to_use <- full_palette[1:length(unique_cell_types)]
      names(palette_to_use) <- unique_cell_types
      
      cytomarker_palette(palette_to_use)
      
      if (isTruthy(markers$top_markers) & length(markers$top_markers) > 1) {
        
        if (input$panel_sorter == "Group by cell type") {
          markers$top_markers <- names(sort(markers$associated_cell_types[names(markers$associated_cell_types) %in%
                                                                            markers$top_markers]))
        } else {
          markers$top_markers <- sort(markers$top_markers)
        }
        
      }
      
      if (isTruthy(markers$scratch_markers) & length(markers$scratch_markers) > 1) {
        if (input$panel_sorter == "Group by cell type") {
          markers$scratch_markers <- names(sort(markers$associated_cell_types[names(markers$associated_cell_types) %in%
                                                                                markers$scratch_markers]))
          
        } else {
          markers$scratch_markers <- sort(markers$scratch_markers)
        }
      }
      
      markers_with_type(markers$associated_cell_types[names(markers$associated_cell_types) %in% 
                                                        markers$top_markers])
      
      # set the marker inclusion list depending on whether or not the configuration
      # is subsetting to Abcam products or not
      # Option 1: if the configuration is true, then have a star for any product in the Abcam catalogue
      # Option 2: if the configuration is false, then have a star if the marker was in the initial suggestions
      marker_filtration <- if(isTruthy(STAR_FOR_ABCAM)) unique(antibody_info$Symbol) else original_panel()
      
      labels_top <- lapply(markers$top_markers, 
                           function(x) div(x, map_gene_name_to_antibody_icon(x, marker_filtration), 
                                           style=paste('padding: 3px; color:', 
                                                       set_text_colour_based_on_background(cytomarker_palette()[ markers$associated_cell_types[x]]), 
                                                       '; background-color:', 
                                                       cytomarker_palette()[ markers$associated_cell_types[x] ])))
      labels_scratch <- lapply(markers$scratch_markers, 
                               function(x) div(x, map_gene_name_to_antibody_icon(x, marker_filtration), 
                                               style=paste('padding: 3px; color:', 
                                                           set_text_colour_based_on_background(cytomarker_palette()[ markers$associated_cell_types[x]]), 
                                                           '; background-color:', 
                                                           cytomarker_palette()[ markers$associated_cell_types[x] ])))
      
      # output$legend <- renderPlot(cowplot::ggdraw(get_legend(cytomarker_palette())))
      
      output$legend <- renderPlot({
        op <- par(family = "sans", mar = c(4, 0.2, 3, 0.2))
        
        longest_factor <- max(nchar(names(cytomarker_palette())))
        label_size <- ifelse(longest_factor < 35, 1.2, ifelse(longest_factor >= 35 & longest_factor < 50,
                                                              1, 0.75))
        plot(NULL ,xaxt='n',yaxt='n',bty='n',
             ylab='',xlab='', xlim=0:1, ylim=0:1)
        legend("top", legend = names(cytomarker_palette()), 
               pch=15, pt.cex=2.2, cex=label_size, bty='n',
               # ncol = ceiling(length(cytomarker_palette())/15),
               ncol = 1,
               col = cytomarker_palette())
        mtext("Cell Type", cex=1.5)
        par(op)
      })
      
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
        
        set_allowed_genes()
        
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
        cells_per_type(table(SummarizedExperiment::colData(sce()[,
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
          
          previous_overall_score(current_overall_score())
        }
        mm <- metrics()
        # if (is.null(mm)) {
        #   return(NULL)
        # }
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
        
        current_overall_score(list(score = round(unique(subset(current_metrics(), what == "Overall")$`score`), 3),
                                   counts = unique(subset(current_metrics(), what == "Overall")$`Counts`)))
        
        current_metrics(list(all = current_metrics() |> filter(what != "Overall"),
                             summary = current_metrics() |>
                               mutate(Counts = as.numeric(Counts)) |>
                               group_by(what, Counts) |>
                               summarize(mean_score = round(mean(score),
                                                            3)) |>
                               arrange(desc(Counts)) |> rename(`Cell Type` = what,
                                                               `Mean Score` = mean_score) |>
                               filter(`Cell Type` != "Overall")))
        
        if (isTruthy(previous_metrics()$all)) {
          all_metrics <- rbind(previous_metrics()$all, current_metrics()$all) |> 
            mutate(what = factor(what,
                                 levels = c(rev(unique(what[what != "Overall"])), "Overall")))
        } else {
          all_metrics <- current_metrics()$all |>  mutate(what = factor(what,
                                                                        levels = c(rev(unique(what[what != "Overall"])), "Overall")))
        }
        
        plots$metric_plot <- suppressWarnings(plotly::plot_ly(all_metrics, x = ~score, y = ~what, 
                                                              color = ~Run,
                                                              type='box', hoverinfo = 'none') %>% 
                                                layout(boxmode = "group",
                                                       xaxis = list(title="Score"),
                                                       yaxis = list(title="Source")))
        
        plots_for_markdown(list(top = plots_for_markdown()$top,
                                all = plots_for_markdown()$all,
                                metric = plots$metric_plot))
        
        previous_run_log_2(previous_run_log())
        previous_run_log(current_run_log())
        
        run_config <- create_run_param_list(marker_list = current_markers(),
                                            input_file = input$input_scrnaseq$datapath,
                                            assay_used = pref_assay(),
                                            het_source = column(),
                                            panel_size = input$panel_size,
                                            cell_cutoff_value = as.integer(cell_min_threshold()),
                                            subsample = input$subsample_sce,
                                            subsample_number = ifelse(isTruthy(input$subsample_sce),
                                                                      input$subset_number, "None"),
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
        
        heatmap_for_report(list(
          correlation = create_heatmap(sce()[,sce()$keep_for_analysis == "Yes"], 
                                       current_markers(), column(), "Marker-marker correlation", 
                                       "Expression", pref_assay()),
          expression = create_heatmap(sce()[,sce()$keep_for_analysis == "Yes"], 
                                      current_markers(), column(), column(), 
                                      "Expression", pref_assay()),
          z_score = create_heatmap(sce()[,sce()$keep_for_analysis == "Yes"], 
                                   current_markers(), column(), column(), 
                                   "z-score", pref_assay())))
        
        
        
        # Show help text popover
        shinyjs::show(id = "marker_visualization")
        shinyjs::show(id = "marker_display")
        
        setProgress(value = 1)
        
        
      })
    }
    
    pre_upload_configuration <- function() {
      updateCheckboxInput(inputId = "precomputed_dim", value = F)
      if (!isTruthy(curated_selection())) use_precomputed_umap(FALSE)
      if (!isTruthy(curated_selection())) umap_precomputed_col(NULL)
      specific_cell_types_selected(NULL)
      
    }
    
    post_upload_configuration <- function(input_sce) {
      input_sce <- detect_assay_and_create_logcounts(input_sce)
      input_sce <- parse_gene_names(input_sce, grch38)
      # input_sce <- remove_confounding_genes(input_sce)
      sce(input_sce)
      
      set_allowed_genes()
      
      shinyjs::hide(id = "download_button")
      
      toggle(id = "analysis_button", condition = isTruthy(sce()))
      # toggle(id = "apps", condition = isTruthy(SUBSET_TO_ABCAM))
      
      input_assays <- c(names(SummarizedExperiment::assays(sce())))
      
      # set the current metrics to NULL
      current_metrics(NULL)
      previous_metrics(NULL)
      current_overall_score(NULL)
      previous_overall_score(NULL)
      
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
                          colnames(SummarizedExperiment::colData(sce()))[!grepl("keep_for_analysis", 
                                                                                colnames(SummarizedExperiment::colData(sce())))][1])
      
      
      updateSelectInput(
        session = session,
        inputId = "coldata_column",
        choices = colnames(SummarizedExperiment::colData(sce()))[!grepl("keep_for_analysis", 
                                                                        colnames(SummarizedExperiment::colData(sce())))],
        selected = selection
      )
      
      specific_cell_types_selected(unique(sce()[[selection]]))
      
      if (!isTruthy(input$coldata_column)) {
        column(colnames(SummarizedExperiment::colData(sce()))[1])
      } else {
        column(input$coldata_column)
      }
      
      aliases_table(NULL)
      gene_aliases_to_show(NULL)
      
      output$table_of_gene_aliases <- DT::renderDataTable(
        gene_aliases_to_show()
      )
      
      
      if (isTruthy(input$bl_top) | isTruthy(current_markers())) {
        proceed_with_analysis(FALSE)
        showModal(reset_option_on_change_modal("uploaded a new dataset"))
        
      }
      
      req(proceed_with_analysis())
      
      possible_umap_dims(detect_umap_dims_in_sce(sce()))
      
      toggle(id = "precomputed", condition = length(possible_umap_dims()) > 0)
      
      if (isTruthy(curated_selection()) & !isTruthy(use_precomputed_umap())) {
        updateCheckboxInput(session, inputId = "precomputed_dim",
                            value = T)
      }
      toggle(id = "gene_alias_tag", condition = isTruthy(gene_aliases_to_show()))
      
    }
    
    set_allowed_genes <- function() {
      
      # TODO: for now, do not use antibody app subsets because of the catalog format
      # if (isTruthy(SUBSET_TO_ABCAM) | length(input$select_aa) > 0) 
      #   allowed_genes(get_allowed_genes(input$select_aa, applications_parsed,
      #   sce()[,sce()$keep_for_analysis == "Yes"])) else
      
      allowed_genes(rownames(sce()[,sce()$keep_for_analysis == "Yes"]))
      
      allowed_genes(if (isTruthy(ONLY_PROTEIN_CODING)) 
        allowed_genes()[allowed_genes() %in% cytomarker_data$protein_coding] else allowed_genes())
      
      allowed_genes(remove_confounding_genes(allowed_genes()))
      
    }
    
    ### SAVE PANEL ###
    output$downloadData <- downloadHandler(
      filename <- paste0("cytomarker-Panel-", Sys.Date(), ".zip"),
      
      content = function(fname) {
        showNotification("Rendering output report and config file, this may take a few moments..", duration = 4)
        # future_promise({
        
        if (!isTruthy(umap_all())) {
          umap_all(get_umap(sce()[,sce()$keep_for_analysis == "Yes"],
                            column(), pref_assay(), 
                            use_precomputed_umap(),
                            umap_precomputed_col(), F, allowed_genes()))
        }
        
        if (!isTruthy(umap_top())) {
          umap_top(get_umap(sce()[,sce()$keep_for_analysis == "Yes"][current_markers()$top_markers,], 
                            column(), pref_assay(), use_precomputed_umap(),
                            umap_precomputed_col(),
                            T, current_markers()$top_markers))
          
        }
        
        plots_for_markdown(list(all = suppressWarnings(plotly::plot_ly(umap_all(), x=~UMAP_1, y=~UMAP_2, color=~get(input$coldata_column), 
                                                                       text=~get(input$coldata_column), 
                                                                       type='scatter', hoverinfo="text", colors=cytomarker_palette()) %>% 
                                                         layout(title = "UMAP all genes")),
                                top = suppressWarnings(plotly::plot_ly(umap_top(), x=~UMAP_1, y=~UMAP_2, color=~get(input$coldata_column), 
                                                                       text=~get(input$coldata_column), 
                                                                       type='scatter', hoverinfo="text", colors=cytomarker_palette()) %>% 
                                                         layout(title = "UMAP selected markers")),
                                metric = plots_for_markdown()$metric))
        
        download_data(fname, current_run_log()$map, plots_for_markdown(), heatmap_for_report(), df_antibody(), 
                      markdown_report_path, current_metrics()$summary, current_overall_score(),
                      markers_with_type())
        # })
      },
      contentType = "application/zip"
    )
    
  }
  
  shinyApp(ui, server, options = list(shiny.autoload.r=FALSE), ...)
}

