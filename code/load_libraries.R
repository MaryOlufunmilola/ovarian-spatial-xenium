# ----------------------------------------------------------
# Loading Libraries
# ----------------------------------------------------------

# ----------------------------------------------------------
# PATHS NOTE FOR LOCAL USE
# This script uses Code Ocean fixed paths (/data/, /results/).
# To run locally, replace with:
#   dataDir    <- here::here("data/")
#   resultsDir <- here::here("results/")
# and add "here" back to cran_pkgs.
# ----------------------------------------------------------

# pacman / BiocManager bootstrap
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# List of CRAN and general-use packages
# NOTE: 'here' removed — it relies on .Rproj which does not exist on Code Ocean
cran_pkgs <- c(
  "devtools", "future", "parallel", "tidyr", "dplyr", "arrow",
  "scCustomize", "SingleR", "Seurat", "ggplot2", "ggtext", "UCell", "ggrepel", "purrr",
  "harmony", "pheatmap", "dbscan", "NMF", "RColorBrewer", "igraph", "svglite",
  "reticulate", "ggnetwork", "CellChat", "presto"
)

# List of Bioconductor-only packages
bioc_pkgs <- c("BiocManager", "BiocParallel", "org.Hs.eg.db")

# Install Bioconductor packages if missing
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = TRUE, ask = FALSE)
}

# Load everything with pacman
pacman::p_load(char = c(cran_pkgs, bioc_pkgs))

# GitHub packages — pre-installed via postInstall; p_load_gh just loads them here
pacman::p_load_gh(c(
  "jokergoo/ComplexHeatmap",
  "jokergoo/circlize",
  "jinworks/CellChat",
  "immunogenomics/presto"
))

# ----------------------------------------------------------
# Paths  (Code Ocean fixed paths — do not change)
# ----------------------------------------------------------
dataDir    <- "/data/"      # upload your Xenium files to /data in the capsule
resultsDir <- "/results/"   # Code Ocean creates /results automatically

# ----------------------------------------------------------
# Options
# ----------------------------------------------------------
options(Seurat.object.assay.version = "v5")
options(
  future.globals.maxSize = 600 * 1024^3,
  future.rng.onMisuse    = "ignore"
)
