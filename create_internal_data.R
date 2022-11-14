
antibody_info <- read_csv(file.path("inst", "Abcam_published_monos_with_gene_v2.csv"))
grch38 <- read_tsv(file.path("inst", "annotables_grch38.tsv"))

all_zones <- OlsonNames()
sl <- grepl("/", all_zones)
all_zones <- all_zones[!sl]

cytosel_data <- list(antibody_info = antibody_info, grch38 = grch38,
                     time_zones = all_zones)

usethis::use_data(cytosel_data, internal = T, overwrite = T)

