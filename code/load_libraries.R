# ----------------------------------------------------------
# Loading Libraries
# ----------------------------------------------------------
# PATHS NOTE FOR LOCAL USE:
# Replace dataDir/resultsDir with here::here() equivalents
# and add "here" back to cran_pkgs.
# ----------------------------------------------------------

# Load all packages, pre-installed at build time, no installation here
pkgs <- c(
  # GitHub
  "ComplexHeatmap", "circlize", "CellChat", "presto",
  # CRAN
  "future", "parallel", "parallelly",
  "tidyr", "dplyr", "arrow",
  "Seurat", "scCustomize",
  "ggplot2", "ggrepel", "purrr",
  "harmony", "pheatmap", "RColorBrewer",
  "igraph", "data.table",
  # Bioconductor
  "BiocParallel", "SingleR", "glmGamPoi", "UCell" 
)

for (pkg in pkgs) {
  library(pkg, character.only = TRUE)
}

# ----------------------------------------------------------
# Paths (Code Ocean fixed paths)
# ----------------------------------------------------------
dataDir    <- "/data/"
resultsDir <- "/results/"

# ----------------------------------------------------------
# Options
# ----------------------------------------------------------
options(Seurat.object.assay.version = "v5")
options(
  future.globals.maxSize = 600 * 1024^3,
  future.rng.onMisuse    = "ignore"
)