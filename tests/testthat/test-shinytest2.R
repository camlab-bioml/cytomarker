library(shinytest2)
library(cytosel)

test_that("{shinytest2} recording: cytosel", {
  testthat::local_edition(3)
  app <- AppDriver$new(cytosel::cytosel(),
                       variant = platform_variant(), name = "cytosel", height = 732, 
      width = 1161, load_timeout = 1e+05
      # shinyOptions = list(test.mode = TRUE)
      )
    announce_snapshot_file("cytosel-001.png")
    skip_if(interactive())
    app$expect_screenshot()
  
})