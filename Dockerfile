FROM rocker/verse
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgeos-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('devtools','BiocManager', 'DT', 'shinyalert', 'reactable', 'bsplus', 'shinyjs', 'plyr', 'shinytest2', 'shinydashboard', 'shinyBS', 'rdrop2', 'shinyanimate'), dependencies=F)"
RUN R -e "BiocManager::install(c('stringr', 'yardstick', 'htmltools', 'methods', 'stats', 'gridExtra', 'testthat', 'tools', 'scran', 'Seurat', 'caret', 'scuttle', 'scater', 'SummarizedExperiment', 'SingleCellExperiment', 'clustifyr', 'zellkonverter', 'printr', 'ggplot2', 'heatmaply', 'plotly', 'geneBasisR', 'waiter', 'cicerone'))"

RUN R -e "devtools::install_github('stephenturner/annotables'); devtools::install_github('camlab-bioml/inferorg'); devtools::install_github('MarioniLab/geneBasisR'); install.packages(c('sortable', 'randomcoloR', emojifont'))"

COPY . cytosel/
RUN cd cytosel

RUN yes | R -e "setwd('cytosel'); options(needs.promptUser = FALSE); library(devtools); devtools::load_all()"  
EXPOSE 3838
CMD ['R', 'devtools::load_all(); cytosel()']
