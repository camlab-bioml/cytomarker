
## Read the tooltips yaml on package load
tooltips <- yaml::read_yaml(
  # system.file("inst", "tooltips.yml", package="cytosel")
  file.path("inst", "tooltips.yml")
)


#' Returns the tooltip as stored
#' in inst/tooltips.yml
get_tooltip <- function(tip_name) {
  tooltips[[tip_name]]
}

