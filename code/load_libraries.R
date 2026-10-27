# ----------------------------------------------------------
# Libraries, paths, and options
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
# Paths: Code Ocean defaults, overridable for local runs
#   DATA_DIR=data/ RESULTS_DIR=results/ CODE_DIR=code Rscript code/analysis.R
# ----------------------------------------------------------
dataDir    <- Sys.getenv("DATA_DIR", "/data/")
resultsDir <- Sys.getenv("RESULTS_DIR", "/results/")
dir.create(resultsDir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# Options
# ----------------------------------------------------------
options(Seurat.object.assay.version = "v5")
options(
  future.globals.maxSize = 600 * 1024^3,
  future.rng.onMisuse    = "ignore"
)