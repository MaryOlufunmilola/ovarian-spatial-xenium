# ----------------------------------------------------------
# Spatial CellChat analysis for the VSTM5 high and low groups
# Authors: Funmi Oyebamiji
# Depends on: Seurat object produced by analysis.R
# Usage:
#   Rscript /code/figure_script_1.R Low
#   Rscript /code/figure_script_1.R High
# ----------------------------------------------------------

# Code Ocean runs from /code; set CODE_DIR=code to run from the repository root.
code_dir <- Sys.getenv("CODE_DIR", "/code")
source(file.path(code_dir, "load_libraries.R"))
source(file.path(code_dir, "functions.R"))

# Read condition from command-line argument ────────────────────────────────
args      <- commandArgs(trailingOnly = TRUE)
condition <- if (length(args) > 0) args[1] else "Low"   # default to Low if not supplied
stopifnot("condition must be 'Low' or 'High'" = condition %in% c("Low", "High"))
message("Running CellChat pipeline for condition: ", condition)

# Load Seurat object ───────────────────────────────────────────────────────
Seurat_object <- readRDS(file.path(resultsDir, "Seurat_object.rds"))
Seurat_object <- UpdateSeuratObject(Seurat_object)
 
# Subset to one condition ──────────────────────────────────────────────────
cells_to_keep <- rownames(Seurat_object@meta.data)[
  Seurat_object@meta.data$VSTM5_Orig_ident == condition
]
seurat_subset <- subset(Seurat_object, cells = cells_to_keep)
 
# Free the full object from RAM before running pipeline ────────────────────
rm(Seurat_object)
gc()
 
# Run CellChat pipeline on single condition ────────────────────────────────
cellchat_result <- run_cellchat_pipeline(
  seurat_list = setNames(list(seurat_subset), condition),
  assay       = "SCT"
)[[condition]]
 
# Free subset from RAM ─────────────────────────────────────────────────────
rm(seurat_subset)
gc()
 
# Save output ──────────────────────────────────────────────────────────────
out_file <- file.path(resultsDir, paste0("vstm5", tolower(condition), ".rds"))
saveRDS(cellchat_result, file = out_file)
message("Saved: ", out_file)