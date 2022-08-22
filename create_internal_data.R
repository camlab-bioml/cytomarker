
antibody_info <- read_csv(file.path("data", "Abcam_published_monos_with_gene.csv"))
grch38 <- read_tsv(file.path("data", "annotables_grch38.tsv"))

cytosel_data <- list(antibody_info = antibody_info, grch38 = grch38)

usethis::use_data(cytosel_data, internal = T, overwrite = T)

