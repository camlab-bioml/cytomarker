library(shinytest2)
library(cytosel)

test_that("{shinytest2} recording: cytosel", {
  app <- AppDriver$new(app_dir = test_path("../../"),
                       variant = platform_variant(), name = "cytosel", height = 732, 
      width = 1161, load_timeout = 1e+05
      # shinyOptions = list(test.mode = TRUE)
      )
  app$expect_screenshot()
})

test_that("{shinytest2} recording: cytosel-select-curated-1", {
  app <- AppDriver$new(app_dir = test_path("../../"),
                       variant = platform_variant(), name = "cytosel-select-curated-1", 
      height = 732, width = 1161, load_timeout = 1e+05)
  app$click("curated_dataset")
  app$click("pick_curated")
  app$set_inputs(curated_options = "PBMC (Blood/Immune)")
  app$expect_screenshot()
})
