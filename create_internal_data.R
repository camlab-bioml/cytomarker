
library(org.Hs.eg.db)
library(annotables)
library(dplyr)
library(feather)
library(tidyverse)
library(stringr)

abcam_antibodies <- read_csv(file.path("inst", "Abcam_published_monos_with_gene_v2.csv"))

grch38 <- read_tsv(file.path("inst", "annotables_grch38.tsv"))

abcam_antibodies <- abcam_antibodies |> dplyr::rename(Symbol = 
                                                        `Gene Name (Upper)`) |>
  tidyr::drop_na()

bd_antibodies <- read_csv(file.path("inst", "bd_us_catalog_20230410_HJackson.csv"))
colnames(bd_antibodies) <- str_to_title(colnames(bd_antibodies))

# alternate method
bd_antibodies <- bd_antibodies |> mutate(Symbol = str_split_fixed(Symbol, " ", 2)[,1]) |>
  mutate(Symbol = str_split_fixed(Symbol, ",|/", 2)[,1])


bd_antibodies <- bd_antibodies |> mutate(Vendor = "BD Biosciences") |> filter(Symbol %in% grch38$symbol) |>
  na.omit()


### View the colnames and modify to fit each other
colnames(abcam_antibodies)
colnames(bd_antibodies)

abcam_antibodies <- abcam_antibodies |> rename(`Product Description` = `Product Category Tier 3`) |>
  mutate(Vendor = "Abcam")
bd_antibodies <- bd_antibodies |> rename(`Datasheet URL` = `Datasheet Url`)

common_cols <- intersect(colnames(bd_antibodies), colnames(abcam_antibodies))

merged_datasets <- merge(abcam_antibodies, bd_antibodies, by=common_cols, all=TRUE)

merged_datasets <- merge(merged_datasets, grch38, 
                         by.x = "Symbol", by.y = "symbol", all.x = T) |>
  select(-c("Ab ID", "Format", "Product Size", "Isotype"))

merged_datasets <- merged_datasets |> relocate(Vendor, .after = `Product Name`)


all_zones <- OlsonNames()
sl <- grepl("/", all_zones)
all_zones <- all_zones[!sl]


applications_parsed <- get_antibody_applications(antibody_info, 
                                                 'Symbol', 'Listed Applications')

gene_mapping <- as.data.frame(unlist(org.Hs.egALIAS2EG))

gene_mapping <- gene_mapping |> group_by(alias_symbol) |> slice_head(n = 1) |> ungroup()

# https://www.genenames.org/download/statistics-and-files/
protein_coding_genes <- (annotables::grch38 |> filter(biotype == "protein_coding"))$symbol

cytomarker_data <- list(# antibody_info = antibody_info,
                      grch38 = grch38,
                     time_zones = all_zones,
                     # applications = applications_parsed,
                     gene_mapping = gene_mapping,
                     protein_coding = protein_coding_genes)

usethis::use_data(cytomarker_data, internal = T, overwrite = T,
                  compress = "xz")

path <- file.path("inst", "merged_catalog.feather")
write_feather(merged_datasets, path)

