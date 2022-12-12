
require(cytosel)

test_that("{shinytest2} recording: cytosel", {
  testthat::local_edition(3)
  
  # set the app depending on the interactive execution to avoid a no package found error
  if (interactive()) {
    cytosel_app <- test_path("../../")
  } else {
    cytosel_app <- cytosel::cytosel()
  }
  app <- AppDriver$new(cytosel_app,
                       variant = platform_variant(), name = "cytosel", height = 732,
      width = 1161, load_timeout = 1e+05
      # shinyOptions = list(test.mode = TRUE)
      )
    announce_snapshot_file("cytosel-001.png")
    # IMPORTANT: only run the tests non-interactively using Github actions
    # currently devtools test and devtools check give different screenshot results
    skip_if(interactive())
    app$expect_screenshot()

})