# cytocel

Antibody panel selection for mass cytometry using scRNA-seq

### Dependencies

To install dependencies:

```r
dependencies <- c(
"scran",
"SingleCellExperiment",
"tibble",
"caret",
"yardstick",
"naivebayes",
"ComplexHeatmap",
"viridis",
"dplyr",
"shiny",
"sortable",
"ggplot2",
"scater",
"forcats",
"cowplot",
"readr"
)
BiocManager::install(dependencies)
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