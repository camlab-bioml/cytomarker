# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-04-20

### Added 

- Config variable `filter_human_gene_names`: yes/no (TRUE/FALSE) to disable filtering for human genes names, enabling 
non-human gene names and cite-seq naming

### Fixed

- Explicitly specify `choices` for certain `selectInput` updates to match new `shiny` API
- Fixed tests for `annotables` gene parsing to be compatible with latest human Ensemble 109

### Removed/Deprecated

- Deprecated: [geneBasisR](https://github.com/MarioniLab/geneBasisR) due to non-active maintenance of the library
and breaking changes in newer versions of R
- Deprecated: screenshot testing using `shinytest22`

## [0.1.0] - 2024-10-09

First public availability