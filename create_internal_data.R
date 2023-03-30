
library(org.Hs.eg.db)
library(annotables)
library(dplyr)
library(feather)

# antibody_info <- read_csv(file.path("inst", "Abcam_published_monos_with_gene_v2.csv"))
antibody_info <- read_tsv(file.path("inst", "registry_with_symbol.tsv"))
grch38 <- read_tsv(file.path("inst", "annotables_grch38.tsv"))

# abcam_antibody <- read_csv(file.path("inst", "Abcam_published_monos_with_gene_v2.csv"))
# 
# abcam_antibody <- abcam_antibody |> dplyr::select(c(`Gene Name (Upper)`, `Listed Applications`)) |>
#   `colnames<-`(c('Symbol', 'Listed Applications'))

# antibody_info <- merge(antibody_info, abcam_antibody, 
#                        by = "Symbol")


all_zones <- OlsonNames()
sl <- grepl("/", all_zones)
all_zones <- all_zones[!sl]


antibody_info <- antibody_info |> tidyr::drop_na()



  # 
  # dplyr::rename(Symbol = 
  #                 `Gene Name (Upper)`) |>

antibody_info <- merge(antibody_info, grch38, 
                       by.x = "Symbol", by.y = "symbol", all.x = T)

# applications_parsed <- get_antibody_applications(antibody_info, 
#                                                  'Symbol', 'Listed Applications')

gene_mapping <- as.data.frame(unlist(org.Hs.egALIAS2EG))

gene_mapping <- gene_mapping |> group_by(alias_symbol) |> slice_head(n = 1) |> ungroup()

# https://www.genenames.org/download/statistics-and-files/
protein_coding_genes <- (annotables::grch38 |> filter(biotype == "protein_coding"))$symbol

cytosel_data <- list(antibody_info = antibody_info, grch38 = grch38,
                     time_zones = all_zones,
                     # applications = applications_parsed,
                     gene_mapping = gene_mapping,
                     protein_coding = protein_coding_genes)

usethis::use_data(cytosel_data, internal = T, overwrite = T,
                  compress = "xz")

path <- file.path("inst", "registry_with_symbol.feather")
write_feather(antibody_info, path)

