# cytosel

![R check](https://github.com/camlab-bioml/cytosel/actions/workflows/check-package.yml/badge.svg)

Antibody panel selection for mass cytometry using scRNA-seq using R Shiny. cytosel is currently hosted
on a [public shinyapps.io server](https://camlab.shinyapps.io/cytosel/)

## Official Documentation

The official documentation for cytosel using Docusaurus can be found [here](https://camlab-bioml.github.io/cytosel-doc/docs/intro)

## Access

cytosel can be accessed through a public shinyapps.io server [here](https://camlab.shinyapps.io/cytosel/)

### Source code access

For those with source code access to cytosel, installation and access can be achieved through the following commands. Ensure that [R](https://cran.r-project.org/) and [RStudio Desktop from Posit](https://posit.co/download/rstudio-desktop/) are installed for your specific operating system. 

```
git clone https://github.com/camlab-bioml/cytosel.git
cd cytosel
R
# inside R console
library(devtools) # or, if devtools isn't installed
install.packages('devtools')
devtools::load_all(); cytosel()
```

This will prompt devtools to install all listed dependencies from the package. 

Additional dependencies in R that are suggested by the package may also need to be explicitly
installed:

```
additional_deps <- c('stringr',
  'yardstick',
  "htmltools',
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
  
  BiocManager::install(dependencies)
```

### Running locally

In the command line, execute the following commands:

```
# navigate to the directory with the cytosel source code
cd cytosel
R -e 'devtools::load_all(); cytosel()'
```

The final line will display a message providing the IP address and local port:

```
Listening on http://127.0.0.1:7665
```

Navigate to that IP address in your browser to run.

Alteratively, in the RStudio Desktop console: 

```
devtools::load_all(); cytosel()
```

### Docker

cytosel can be run in a containerized environment using the following commands. Ensure that Docker is installed for the specific local OS:

```
cd cytosel
docker build -t cytosel . 
docker run -p 3888:3888 docker.io/library/cytosel 
````

This should prompt the user to navigate to a local URL such as `http://localhost:3888/`


## Developer guide

The developer guide for cytosel is in progress and being hosted [here](https://github.com/camlab-bioml/cytosel/wiki/cytosel-Developer-guide)

