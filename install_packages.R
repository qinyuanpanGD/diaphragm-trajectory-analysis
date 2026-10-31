#!/usr/bin/env Rscript
packages <- c("dplyr", "ggplot2", "readr", "tidyr", "pheatmap", "patchwork")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
