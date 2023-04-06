
curated_datasets <- readr::read_csv(system.file("ts_datasets.csv", package = "cytomarker"))

curated_datasets$tissue

compartments <- list(`Endothelial` = curated_datasets$tissue[! curated_datasets$tissue %in% 
                                                             c("Blood", "Bone_Marrow")],
                     `Epithelial` = curated_datasets$tissue[! curated_datasets$tissue %in% 
                                                               c("Blood", "Bone_Marrow",
                                                                 "Fat", "Lymph_Node",
                                                                 "Muscle", "Spleen")],
                     `Immune` = curated_datasets$tissue,
                     `Stromal` = curated_datasets$tissue[! curated_datasets$tissue %in% 
                                                               c("Blood", "Bone_Marrow",
                                                                 "Kidney", "Spleen")]
                     
                     )

write_yaml(compartments, file.path("inst", "ts_compartments.yml"))
