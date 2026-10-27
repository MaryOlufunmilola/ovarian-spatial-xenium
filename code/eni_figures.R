# ----------------------------------------------------------
# QC and Diagnostic Figures
# Author: Funmi Oyebamiji
# Depends on: Seurat_object.rds from analysis.R
# Diagnostic figures, not in the paper (paper figures come from figure_script_2.R).
# Optional: not called by `run`. Rscript code/eni_figures.R
# ----------------------------------------------------------

# Code Ocean runs from /code; set CODE_DIR=code to run from the repository root.
code_dir <- Sys.getenv("CODE_DIR", "/code")
if (!exists("dataDir"))      source(file.path(code_dir, "load_libraries.R"))
if (!exists("save_dimplot")) source(file.path(code_dir, "functions.R"))

# Load objects saved by analysis.R
xen_integrated <- readRDS(file.path(resultsDir, "Seurat_object.rds"))

# Violin plots for each sample before filtering are written by analysis.R (Step 1).

# ── Unintegrated UMAP — coloured by sample ───────────────────────────────────
save_dimplot(xen_integrated,
             filename  = file.path(resultsDir, "DimPlot.umap.unintegrated.png"),
             reduction = "umap.unintegrated", group_by = "sample",
             cols = sample_colors, raster = FALSE, label = TRUE)

# ── Harmony UMAP — coloured by cluster ───────────────────────────────────────
save_dimplot(xen_integrated,
             filename  = file.path(resultsDir, "DimPlot.umap.original.png"),
             reduction = "umap", group_by = "runharmony_clusters",
             cols = cluster_colors, raster = FALSE, label = TRUE)

# ── Subclustered UMAP ────────────────────────────────────────────────────────
save_dimplot(xen_integrated,
             filename  = file.path(resultsDir, "DimPlot.umap.combined.png"),
             reduction = "umap", group_by = "runharmony_combined",
             cols = cluster_colors, raster = FALSE, label = TRUE)

# ── SingleR annotation UMAPs ─────────────────────────────────────────────────
for (ref_name in c("hpca", "blueE")) {
  save_dimplot(xen_integrated,
               filename   = file.path(resultsDir, paste0("DimPlot.umap.", ref_name, ".png")),
               reduction  = "umap",
               group_by   = paste0(ref_name, "_cellmain"),
               cols       = singler_colors,
               raster     = FALSE,
               label      = TRUE,
               label_size = 3)
}

# ── Cell type UMAP — with legend (only version used) ─────────────────────────
umap_theme <- list(
  ggplot2::ggtitle("UMAP by Cell Type"),
  ggplot2::xlab("UMAP1"), ggplot2::ylab("UMAP2"),
  ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5))
)
p <- DimPlot(xen_integrated, reduction = "umap", group.by = "Celltypes",
             cols = cell_colors, raster = FALSE, label = TRUE, label.size = 4) +
     umap_theme
ggplot2::ggsave(file.path(resultsDir, "DimPlot.umap.celltypes.WiLegend.png"),
                p, height = 9, width = 13, units = "in", dpi = 300)

message("All QC figures saved to ", resultsDir)
