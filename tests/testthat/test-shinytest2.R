require(cytomarker)

test_that("{shinytest2} recording: cytomarker", {
  testthat::local_edition(3)

  # set the app depending on the interactive execution to avoid a no package found error
  if (interactive()) {
    cytomarker_app <- test_path("../../")
  } else {
    cytomarker_app <- cytomarker::cytomarker()
  }
  # cytomarker_app <- test_path("../../")
  announce_snapshot_file("cytomarker-001.png")
  # announce_snapshot_file("cytomarker-002.png")
  app <- AppDriver$new(cytomarker_app,
                       variant = platform_variant(), name = "cytomarker", height = 732,
      width = 1161, load_timeout = 1e+05
      # shinyOptions = list(test.mode = TRUE)
      )
    
    # Sys.sleep(7)
    app$wait_for_value(output = "cytomarker_logo")
    # IMPORTANT: only run the tests non-interactively using Github actions
    # currently devtools test and devtools check give different screenshot results
    
    skip_if(interactive())
    
    app$expect_screenshot()
    
    # skip_if(interactive())
    # app$click("help_guide")
    # Sys.sleep(7)
    # app$expect_screenshot()
    
})
