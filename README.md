# cytosel

![R check](https://github.com/camlab-bioml/cytosel/actions/workflows/check-package.yml/badge.svg)

Antibody panel selection for mass cytometry using scRNA-seq

### To run

Navigate to the root directory and run

```r
devtools::load_all(); cytosel()
```

### Dependencies

To install dependencies:

```r
dependencies <- c(
"scran",
"SingleCellExperiment",
"tibble",
"caret",
"naivebayes",
"dplyr",
"shiny",
"sortable",
"ggplot2",
"scater",
"forcats",
"cowplot",
"readr",
"dqshiny",
"DT",
"shinyalert",
"SummarizedExperiment",
"scuttle",
"Seurat",
"viridis",
"reactable",
"ComplexHeatmap",
"bsplus",
"shinyjs",
"zip"
)
BiocManager::install(dependencies)
remotes::install_github("daqana/dqshiny")
```

### To run

In the command line, run

```bash
git clone https://github.com/camlab-bioml/cytosel.git
R -e 'shiny::runApp("cytosel")'
```

The final line will display something like

```
Listening on http://127.0.0.1:7665
```

Navigate to that IP address in your browser to run.
