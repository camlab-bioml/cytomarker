library(shinytest2)

app <- AppDriver$new(test_path("../../"), load_timeout = 100000)
record_test(app, load_timeout = 100000)
