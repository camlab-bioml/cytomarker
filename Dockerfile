FROM rocker/verse
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgsl-dev \
    libgeos-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('devtools','BiocManager', 'DT', 'shinyalert', 'reactable', 'bsplus', 'shinyjs', 'plyr', 'shinytest2', 'shinydashboard', 'shinyBS', 'rdrop2', 'shinyanimate', 'stringr', 'yardstick', 'htmltools', 'methods', 'stats', 'gridExtra', 'testthat', 'tools', 'printr', 'ggplot2', 'heatmaply', 'plotly', 'waiter', 'cicerone'), dependencies=T)"
RUN R -e "BiocManager::install(c('scran', 'Seurat', 'caret', 'scuttle', 'scater', 'SummarizedExperiment', 'SingleCellExperiment', 'clustifyr', 'zellkonverter'))"


RUN R -e "devtools::install_github('stephenturner/annotables', dependencies = F); devtools::install_github('camlab-bioml/inferorg', dependencies = F)"

RUn R -e "install.packages(c('sortable', 'randomcoloR', 'emojifont'))"

RUN R -e "devtools::install_github('MarioniLab/geneBasisR', dependencies = T)"

RUN R -e "devtools::install_github('react-R/reactR')"

RUN mkdir /cytosel

COPY . /cytosel

COPY R/ /cytosel/R

WORKDIR /cytosel

RUN ls

# RUN cd cytosel

# RUN yes | R -e "setwd('cytosel'); options(needs.promptUser = FALSE); library(devtools); devtools::load_all()"  


ENV PORT=3838

CMD R -e "options(shiny.port = 3838, shiny.host='0.0.0.0'); devtools::load_all(); cytosel()"
