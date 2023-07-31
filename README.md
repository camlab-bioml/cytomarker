# cytomarker

![R check](https://github.com/camlab-bioml/cytomarker/actions/workflows/check-package.yml/badge.svg)

Antibody panel selection for mass cytometry using scRNA-seq using R Shiny. cytomarker is currently hosted
on a [public shinyapps.io server](https://camlab.shinyapps.io/cytomarker/)

## Official Documentation

The official documentation for cytomarker using Docusaurus can be found [here](https://camlab-bioml.github.io/cytomarker-doc/docs/intro)

## Access

cytomarker can be accessed through a public shinyapps.io server [here](https://camlab.shinyapps.io/cytomarker/)


### Source code access

For those with source code access to cytomarker, installation and access can be achieved through the following commands. Ensure that [R](https://cran.r-project.org/) and [RStudio Desktop from Posit](https://posit.co/download/rstudio-desktop/) are installed for your specific operating system. 

```
git clone https://github.com/camlab-bioml/cytomarker.git
cd cytomarker
```

Once cloned, open RStudio and import the project using `File -> Open Project`, and select the `cytomarker.Rproj` file. After importing the codebase, install the app dependencies:

```
library(devtools) # or, if devtools isn't installed
install.packages('devtools')
devtools::load_all(); cytomarker()
```

This will prompt devtools to install all listed dependencies from the package. 

Additional dependencies in R that are suggested by the package may also need to be explicitly
installed:

```
additional_deps <- c('stringr',
  'yardstick',
  'htmltools',
  'methods',
  'stats',
  'gridExtra',
  'testthat',
  'tools',
  'scran',
  'Seurat',
  'caret',
  'scuttle',
  'scater',
  'SummarizedExperiment',
  'SingleCellExperiment',
  'clustifyr',
  'zellkonverter',
  'printr',
  'ggplot2',
  'heatmaply',
  'plotly')
  
  BiocManager::install(additional_deps)
```

### Running locally

In the command line, execute the following commands:

```
# navigate to the directory with the cytomarker source code
cd cytomarker
R -e 'devtools::load_all(); cytomarker()'
```

The final line will display a message providing the IP address and local port:

```
Listening on http://127.0.0.1:7665
```

Navigate to that IP address in your browser to run.

Alteratively, in the RStudio Desktop console: 

```
devtools::load_all(); cytomarker()
```

## Developer guide

The developer guide for cytomarker is in progress and being hosted [here](https://github.com/camlab-bioml/cytomarker/wiki/cytomarker-Developer-guide)

