#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This sets up the analysis section for the F1000 paper.

## ---- init-analysis ----
.cran_packages  <-  c("knitr", "phyloseqGraphTest", "phyloseq", "shiny",
                    "miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest",
                    "vegan", "plyr", "dplyr", "ggrepel", "nlme",
                    "reshape2", "devtools", "PMA","structSSI","ade4",
                    "igraph", "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("phyloseq", "genefilter", "impute")

# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

.inst <- .github_packages %in% installed.packages()
if (any(!.inst)) {
  devtools::install_github(.github_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}

## ---- ggplot-theme ----
library("ggplot2")
theme_set(theme_bw())
min_theme <- theme_update(panel.border = element_blank(),
                          panel.grid = element_blank(),
                          axis.ticks = element_blank(),
                          legend.title = element_text(size = 8),
                          legend.text = element_text(size = 6),
                          axis.text = element_text(size = 6),
                          axis.title = element_text(size = 8),
                          strip.background = element_blank(),
                          strip.text = element_text(size = 8),
                          legend.key = element_blank())

## ---- misc-setup ----
set.seed(10)
options(width=100)

## ---- get-raw-data ----
raw_data <- function() {
  ps_file <- "data/ps.rds"
  readRDS(ps_file)
}

## ---- get-preprocessed-data ----
preprocessed_data <- function() {
  ps <- raw_data()
  ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
  sample_data(ps)$librarysize <- log10(rowSums(otu_table(ps)))
  sample_data(ps)$age_binned <- cut(sample_data(ps)$age, breaks = c(0, 100, 200, 400))
  ps0 <- ps
  pslog <- transform_sample_counts(ps, function(x) log(1 + x))
  # detected in preprocessing.R
  outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
  ps0 <- prune_samples(!(sample_names(ps0) %in% outliers), ps0)
  pslog <- prune_samples(!(sample_names(pslog) %in% outliers), pslog)
  list(ps0 = ps0, ps = ps, pslog = pslog)
}

## ---- setup-example ----
setup_example <- function(pkgs) {
  # clear workspace, except data getters
  all_obj <- ls(envir = .GlobalEnv)
  all_obj <- setdiff(all_obj, c("raw_data", "preprocessed_data", "setup_example"))
  rm(list = all_obj, envir = .GlobalEnv)
  # add packages and data to workspace
  sapply(pkgs, require, character = TRUE)
  # slight problem with phylo object class...
  attach(preprocessed_data()) # from data.R
}
