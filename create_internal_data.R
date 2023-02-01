
antibody_info <- read_csv(file.path("inst", "Abcam_published_monos_with_gene_v2.csv"))
grch38 <- read_tsv(file.path("inst", "annotables_grch38.tsv"))

all_zones <- OlsonNames()
sl <- grepl("/", all_zones)
all_zones <- all_zones[!sl]


antibody_info <- antibody_info |> dplyr::rename(Symbol = 
                  `Gene Name (Upper)`) |>
  tidyr::drop_na()

antibody_info <- merge(antibody_info, grch38, 
                       by.x = "Symbol", by.y = "symbol", all.x = T)

applications_parsed <- get_antibody_applications(antibody_info, 
                                                 'Symbol', 'Listed Applications')

cytosel_data <- list(antibody_info = antibody_info, grch38 = grch38,
                     time_zones = all_zones, applications = applications_parsed)

usethis::use_data(cytosel_data, internal = T, overwrite = T,
                  compress = "xz")

