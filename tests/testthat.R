if (!interactive()) {
  Sys.setenv("R_TESTS" = "")
}

# Sys.setenv("R_TESTS_VERBOSE" = "2")
# Sys.setenv("R_MAX_MEM_SIZE" = "400M")

library(testthat)
library(shinytest2)
library(cytosel)
library(yaml)
library(S4Vectors)

test_check("cytosel")