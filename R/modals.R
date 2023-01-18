#' Show a shinyalert error if the input file type is invalid
#' @importFrom shinyalert shinyalert
invalid_modal <- function() { # Input file is invalid
  shinyalert(
    title = "Error",
    text = paste("Input must be a Single Cell Experiment or Seurat object. Please upload a different file."),
    type = "error",
    showConfirmButton = TRUE,
    confirmButtonCol = "#337AB7"
  )
}

#' #' Show an input modal to select the single cell assay to use for analysis
#' #' @importFrom shiny modalDialog
#' assay_modal <- function(assays, failed = FALSE) { # Assay selection
#'   modalDialog(
#'     selectInput("assay",
#'                 "Choose which assay to load",
#'                 assays),
#'     helpText("Recommended assay type is logcounts, as otherwise panel selection 
#'                  will be skewed towards high abundance transcripts rather than heterogeneously expressed transcripts."),
#'     if (failed) {
#'       div(tags$b("Error", style = "color: red;"))
#'     },
#'     footer = tagList(
#'       actionButton("assay_select", "Select")
#'     )
#'   )
#' }
#' 
#' 


#' Show a shinyalert error if either 1. the selected column has only one type, or
#' 2. the selected column has over 100 unique types
#' @importFrom shinyalert shinyalert
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

#' Show a shinyalert warning if the uploaded marker is not found in the sce used for analysis
#' @importFrom shinyalert shinyalert
warning_modal <- function(not_sce) { # Uploaded marker is not in SCE
  shinyalert(title = "Warning", 
             text = paste(paste(not_sce, collapse = ", "), " are not in SCE."), 
             type = "warning", 
             showConfirmButton = TRUE,
             confirmButtonCol = "#337AB7")
}

#' Show an input modal to select the markers to remove from the selected marker space to scratch
#' @importFrom shiny modalDialog
suggestion_modal <- function(failed = FALSE, suggestions, possible_removal) { # Marker removal suggestion
  modalDialog(
    selectInput("markers_to_remove",
                "Select markers to remove",
                choices = possible_removal,
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

#' Show a shinyalert error if the metadata column selected has too many levels
#' @importFrom shinyalert shinyalert
#' @param column the string name of the metadata column selected
invalid_metadata_modal <- function(column) {
  shinyalert(title = "Error", 
             text = HTML(paste("Metadata column:",'<br/>',
                          "<b>", toString(column), "</b>", '<br/>', " has more than 100 unique levels. Please select another metadata column to analyze.")), 
             type = "error", 
             showConfirmButton = TRUE,
             confirmButtonCol = "#337AB7",
             html = T)
}

#' Show a shinyalert error if the input marker is not found in the current antibody set
#' @importFrom shinyalert shinyalert
#' @param dne the string name of the selected marker
dne_modal <- function(dne) { # Input marker is not in the dataset
  shinyalert(title = "Error", 
             text = paste("Marker ", paste(dne), " does not exist or has no corresponding antibody."), 
             type = "error", 
             showConfirmButton = TRUE,
             confirmButtonCol = "#337AB7")
}

#' Show an input modal to view the frequency tally and proportion count for a selected heterogeneity category.
#' Furthermore, specify a custom subset of cell types to retain for analysis
#' @importFrom shiny modalDialog
#' @param cell_cat_value The minimum cell category cutoff (minimum should be set to 2)
#' @param cell_choices A vector of all possible cell types in the selected category
#' @param cell_types_included vector of cell types that are selected from the category
cell_cat_modal <- function(cell_cat_value, cell_choices, cell_types_included) {
  modalDialog(
    DT::dataTableOutput("cell_cat_table"),
    flowLayout(cellArgs = list(style = c("width: 300px;")),
               selectInput("user_selected_cells", 
                           "Create a custom subset for analysis", 
                           choices = cell_choices, 
                           selected = cell_types_included,
                           multiple=TRUE),
               numericInput("min_category_count", "Minimum cell category cutoff:", 
                            cell_cat_value, min = 2, max = 100, step = 0.5, width = NULL)),
    title = "Frequency Count for selected heterogeneity category",
    helpText("Certain cell types can be manually ignored by the user during analysis in the first dialog box above. If this cell is left empty, then by default 
                 all cell types with a Freq of 2 or greater will be retained for analysis. Additionally, 
                 cell categories below a minimum count threshold can be excluded."),
    size = "xl",
    easyClose = TRUE,        
    footer = tagList(
      actionButton("add_selected_to_analysis", "Set subset for analysis."),
      modalButton("Cancel.")
    ))
}


#' Show a shinyalert error if the minimum cell cutoff input is set below 2
#' @importFrom shinyalert shinyalert
threshold_too_low_modal <- function() { # Input marker is not in the dataset
  shinyalert(title = "Error", 
             text = "Cytosel requires the minimum cell type category to be set to 2 or greater for statistical inference. Please adjust the value of minimum cell category cutoff to at least 2.", 
             type = "error", 
             showConfirmButton = TRUE,
             confirmButtonCol = "#337AB7")
}

#' Show a shinyalert warning if cell types retained for analysis are below the minimum cell cutoff. 
#' cytosel will ignore these for analysis and print their names during analysis.
#' @importFrom shinyalert shinyalert
#' @param min_cat_count The current value of minimum cell cutoff (minimum of 2).
#' @param the vector of the names of cell types below the cutoff to be ignored.
cell_type_ignored_modal <- function(min_cat_count, excluded_cells) {
  shinyalert(title = "Warning",
             text = HTML(paste("Cell types with abundance below the set threshold of ",
                               min_cat_count,
                               ":", '<br/>',
                               "<b>", toString(excluded_cells), "</b>", '<br/>',
                               "were identified and not removed by the user.",
                               "They are to be ignored during analysis. The threshold can be changed using Minimum cell category cutoff.")),
             type = "warning",
             showConfirmButton = TRUE,
             confirmButtonCol = "#337AB7",
             html = TRUE)
}

#' Show a shinyalert error if, after user subsetting and filtering, no cells remain for analysis
#' @importFrom shinyalert shinyalert
no_cells_left_modal <- function() { # Uploaded marker is not in SCE
  shinyalert(title = "Error", 
             text = "No cells remain after filtering/subsetting. Please review the input parameters.", 
             type = "error", 
             showConfirmButton = TRUE,
             confirmButtonCol = "#337AB7")
}

#' Upload markers from a .txt file, add various markers by name, or view suggested markers from a cell category.
#' @param markers_addable The vector of gene markers that can be added to the input. Should be generated using alowed_genes() in the app
#' @param suggest_cell_types The vector of cell types for marker suggestion. Should match the names of the findMarkers() output. 
#' @param session The shiny session currently being used
markers_add_modal <- function(markers_addable, suggest_cell_types, session) { 
  updateSelectizeInput(
    session = session,
    inputId = "add_markers",
    # choices = current_markers()$top_markers,
    # selected = current_markers()$top_markers[1])
    choices = markers_addable,
    selected =markers_addable[1],
    server = T)
  
  modalDialog(
    helpText("Upload markers from a .txt file, add various markers by name, or view suggested markers from a cell category."),
    fileInput("uploadMarkers", "Upload markers", width = "100%"),
    div(style = "margin-bottom: -30px"),
    flowLayout(cellArgs = list(
      style = "margin-top: 0px;
                         margin-right: 0px;
                         margin-bottom: 0px; 
          margin-left: 0px; "), actionButton("add_to_selected", "Add uploaded", width = "100%"),
      actionButton("replace_selected", "Replace selected", width = "100%")),
    br(),
    br(),
    flowLayout(cellArgs = list(
      style = "margin-top: 0px;
                         margin-right: 0px;
                         margin-bottom: -15px; 
          margin-left: 0px; "), 
      selectizeInput("add_markers", "Add marker by name", choices = NULL, width = "100%", multiple = T),
      selectInput('cell_type_markers', "Suggest markers for cell type:", choices=suggest_cell_types)),
    # div(style="display:inline-block",selectizeInput("add_markers", "Add marker by name", 
    #       choices = NULL, width = "100%", multiple = F)),
    # div(style="display:inline-block", selectInput('cell_type_markers', "Suggest markers for cell type:", choices=NULL))
    flowLayout(actionButton("enter_marker", "Add"),
               actionButton('add_cell_type_markers', "Suggest")))
  # div(style="display:inline-block", actionButton("enter_marker", "Add")),
  # div(style="display:inline-block", actionButton('add_cell_type_markers', "Suggest")))
}


#' Show a shinyalert error if the user tries to render the distribution table with more than 100 unique elements
#' @importFrom shinyalert shinyalert
too_large_to_show_modal <- function(col_name) {
  shinyalert(
    title = "Error",
    text = paste("Column", col_name, "has more than 100 unique elements. Please select another column."),
    type = "error",
    showConfirmButton = TRUE,
    confirmButtonCol = "#337AB7"
  )
}


#' Show an input modal if the user wants to use precomputed UMAP coordinates
#' @importFrom shiny modalDialog
pre_computed_umap_modal <- function(possible_dims) { # Marker removal suggestion
  modalDialog(
    title = "Use Precomputed UMAP",
    helpText("Select a possible UMAP dimension detected in the uploaded SCE."),
    flowLayout(
    selectInput("possible_precomputed_dims",
                NULL,
                choices = possible_dims,
                selected = NULL,
                multiple = F),
    actionButton("select_precomputed_umap", "Use selected for UMAP")
    )
    )
}


#' Show a shinyalert warning if the sce for analysis has more than a set threshold of cell counts
#' @importFrom shinyalert shinyalert
#' @param ncell_sce The number of cells in the Singlecellexperiment object used for analysis
#' @param num_cells The threshold for cell number for the sce. An ace with a cell number greater than this threshold will trigger the warning.
high_cell_number_warning <- function(ncell_sce, num_cells) { # Uploaded marker is not in SCE
  
  if (ncell_sce > num_cells) {
    shinyalert(title = paste("Warning: cell counts greater than ", num_cells, sep = ""), 
               text = paste("The object retained for analysis has", 
                            ncell_sce, " cells which is greater than the recommended threshold of:",
                            num_cells, "Consider using the downsample function for better performance.", 
                            sep = " "), 
               type = "warning", 
               showConfirmButton = TRUE,
               confirmButtonCol = "#337AB7")
  }
}

#' Show a shinyalert warning if the dimensions of the current sce do not match the input yaml
#' @importFrom shinyalert shinyalert
#' @param title_message The type of warning message in the title
#' @param body_message The type of warning message in the body
reupload_warning_modal <- function(title_message, body_message) {
  shinyalert(title = paste("Warning: ", title_message, sep = ""), 
             text = HTML(paste("The following parameters in the uploaded yml are different from the current SCE. They will be disregarded: ",
                               '<br>',
                               "<b>", 
                          body_message,
                          "</b>")), 
             type = "warning", 
             showConfirmButton = TRUE,
             confirmButtonCol = "#337AB7",
             html = T)
}

#' Show an input modal to confirm resetting the analysis
#' @importFrom shiny modalDialog
reset_analysis_modal <- function() {
  modalDialog(helpText("Resetting the marker panel will erase the current top marker and scratch marker lists, and allow a new set to be generated. Please confirm that you would like to override the current panel."),
  actionButton("reset_marker_panel", "Confirm reset marker panel"),
  footer = tagList(
    modalButton("Cancel")))
}

#' Show an input modal if the current panel contains genes that are not valid (either not in the dataset or not permitted genes)
#' @importFrom shinyalert shinyalert
#' @param missing_genes vector of genes that are not in the current SCE but in the current marker panel
#' @param location Where the 
current_pan_not_valid_modal <- function(missing_genes, location) { # Marker removal suggestion
  shinyalert(title = "Error",
             text = HTML(paste("The current panel contains the following genes that are not found in the ", location,
                               '<br>',
                               "<b>", toString(missing_genes), "</b>", '<br/>',
                               "Please reset the marker panel.")),
             type = "error",
             showConfirmButton = TRUE,
             confirmButtonCol = "#337AB7",
             html = TRUE)
  
}


#' Show an input modal for the user to select a pre-curated cytosel dataset
#' @importFrom shiny modalDialog
#' @param compartment_options a vector of the possible cell type supersets as established by Tabula Sapiens
#' @param compartments_selected a vector of currently selected cell type supersets
#' @param dataset_options a vector of the identifiers for the possible loadable datasets
curated_dataset_modal <- function(compartment_options, compartments_selected,
                                  dataset_options, failed = FALSE) {
  modalDialog(
    flowLayout(selectInput("curated_options",
                       "1. Choose a tissue type to analyze",
                       dataset_options),
               selectInput("curated_compartments",
                          "2. Select cell compartment(s)",
                          compartment_options, multiple = T, 
                          selected = compartments_selected),
    htmlOutput("curated_set_preview")),
    if (failed) {
      div(tags$b("Error", style = "color: red;"))
    },
    footer = tagList(
      actionButton("pick_curated", "Select"),
      modalButton("Cancel")
    )
  )
}

#' Show an error modal if, after subsampling, certain cell types do not meet the minimum count threshold
#' @importFrom shinyalert shinyalert
#' @param cell_types A vector of the cell types that do not pass the count threshold
subsampling_error_modal <- function(cell_types) {
  shinyalert(title = "Error", 
             HTML(paste("After subsampling, the following cell types have either counts below 2 or proportions that are less than 0.5% of the dataset",
                        ":", '<br/>',
                        "<b>", toString(cell_types), "</b>", '<br/>',
                        "Please filter these cell types from the analysis and re-run.")),
             type = "error", 
             showConfirmButton = TRUE,
             confirmButtonCol = "#337AB7",
             html = TRUE)
}

#' Allow the user to change the time zone in the modal
#' @importFrom shiny modalDialog
time_zone_modal <- function(possible_time_zones, current_input) {
  
  time_zone <- ifelse(isTruthy(current_input), current_input, "EST")
  modalDialog(
    flowLayout(cellArgs = list(
      style = "margin-top:10x;
                         margin-right: 15px;
                         margin-bottom: 0px; 
          margin-left: 15px; "), selectInput("time_zone_options",
                                             "Select a time zone to set",
                                             possible_time_zones,
                                             selected = time_zone),
      htmlOutput("current_time", style = "width: 140%; margin-top: 20px; margin-left:15px;")
  ),
  footer = tagList(
    actionButton("pick_time_zone", "Set time zone"),
    modalButton("Cancel")
  ))
}

#' Show an input modal to give the user the option to reset the panel on 
#' @importFrom shiny modalDialog
#' @param change The type of change detected in the configuration
reset_option_on_change_modal <- function(change) {
  modalDialog(helpText(HTML(paste("You have ", "<b>", change, "</b>", ". The current panel may be invalid with these data. 
                                  Would you like to reset the current panel?", "<br/>"))),
              actionButton("reset_marker_panel_reupload", "Reset the marker panel"),
              footer = tagList(
                actionButton("dismiss_marker_reset", "Dismiss")))
}


