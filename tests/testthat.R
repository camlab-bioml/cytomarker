if (!interactive()) {
  Sys.setenv("R_TESTS" = "")
}

library(testthat)
library(shinytest2)
library(cytosel)
library(yaml)

test_check("cytosel")

shinytest2::test_app()
