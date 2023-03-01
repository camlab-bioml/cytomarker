require(cytosel)

test_that("{shinytest2} recording: cytosel", {
  testthat::local_edition(3)

  # set the app depending on the interactive execution to avoid a no package found error
  if (interactive()) {
    cytosel_app <- test_path("../../")
  } else {
    cytosel_app <- cytosel::cytosel()
  }
  
  announce_snapshot_file("cytosel-001.png")
  announce_snapshot_file("cytosel-002.png")
  app <- AppDriver$new(cytosel_app,
                       variant = platform_variant(), name = "cytosel", height = 732,
      width = 1161, load_timeout = 1e+05
      # shinyOptions = list(test.mode = TRUE)
      )
    
    # Sys.sleep(7)
    app$wait_for_value(output = "cytosel_logo")
    # IMPORTANT: only run the tests non-interactively using Github actions
    # currently devtools test and devtools check give different screenshot results
    skip_if(interactive())
    app$expect_screenshot()

    skip_if(interactive())
    app$click("help_guide")
    Sys.sleep(7)
    app$expect_screenshot()

})
