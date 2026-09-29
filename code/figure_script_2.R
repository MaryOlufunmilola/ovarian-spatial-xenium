# ----------------------------------------------------------
# Dot plots, bar charts, UMAP, CellChat heatmaps
# Figures using Seurat_object.rds and CellChat objects
# Written by Eleanor Paskus, modified by Funmi Oyebamiji
# ----------------------------------------------------------

# Load libraries first
source("/code/load_libraries.R")
source("/code/functions.R")

library(ggtext)
library(officer)
library(rvg)
library(tidyverse)

# ==========================================================================
# Figure utility functions
# ==========================================================================

#' Save a ggplot or ComplexHeatmap object to a PowerPoint file
#'
#' @param plot_obj   A ggplot object, or a ComplexHeatmap (set is_complex = TRUE)
#' @param filename   Output filename only e.g. "myplot.pptx", saved to /results/
#' @param width      Placeholder width in inches
#' @param height     Placeholder height in inches
#' @param left       Left offset in inches (default 0)
#' @param top        Top offset in inches (default 0)
#' @param is_complex TRUE for ComplexHeatmap objects that require draw()
#' @param slide_w    Retained for compatibility, slide resizing done manually in PowerPoint
#' @param slide_h    Retained for compatibility, slide resizing done manually in PowerPoint
save_pptx <- function(plot_obj, filename, width, height,
                      left = 0, top = 0, is_complex = FALSE,
                      slide_w = NULL, slide_h = NULL) {
  doc <- officer::read_pptx()
  doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
  value <- if (is_complex) {
    rvg::dml(code = { ComplexHeatmap::draw(plot_obj) })
  } else {
    rvg::dml(ggobj = plot_obj)
  }
  doc <- officer::ph_with(
    x        = doc,
    value    = value,
    location = officer::ph_location(left = left, top = top, width = width, height = height)
  )
  print(doc, target = file.path(resultsDir, filename))
  message("Saved: ", file.path(resultsDir, filename))
}

#' Viridis-filled DotPlot with outlined circles
#'
#' @param obj          Seurat object
#' @param features     Character vector of gene names
#' @param fill_limits  Numeric vector c(min, max) for fill scale
#' @param scale        Passed to DotPlot scale argument (default FALSE)
styled_dotplot <- function(obj, features, fill_limits, scale = FALSE) {
  p <- Seurat::DotPlot(object = obj, features = features, scale = scale) +
    ggplot2::scale_size(range = c(1, 8)) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = NULL, y = NULL)
  p$layers[[1]]$mapping$fill      <- p$layers[[1]]$mapping$colour
  p$layers[[1]]$mapping$colour    <- NULL
  p$layers[[1]]$aes_params$colour <- "black"
  p$layers[[1]]$aes_params$shape  <- 21
  p$layers[[1]]$aes_params$stroke <- 0.5
  p + ggplot2::scale_fill_viridis_c(
    limits = fill_limits,
    name   = "Average Expression",
    oob    = scales::squish
  ) + ggplot2::guides(color = "none")
}


# ==========================================================================
# Load objects
# ==========================================================================

# Load objects saved by earlier scripts
Seurat_object <- readRDS(file.path(resultsDir, "Seurat_object.rds"))
vstm5low      <- readRDS(file.path(resultsDir, "vstm5low.rds"))
vstm5high     <- readRDS(file.path(resultsDir, "vstm5high.rds"))

# Add new_ifit column from IFIT1B_Orig_ident
Seurat_object@meta.data[["new_ifit"]] <- Seurat_object@meta.data[["IFIT1B_Orig_ident"]]

# ==========================================================================
# SECTION 1, Plots using Seurat_object
# ==========================================================================

# VSTM5 & IFIT1B expression per cluster 
# Manual steps automated:
#   - legend order: Percent Expressed first, Average Expression second
#   - shaded percent expressed circles in legend
#   - bold legend category labels (ggtext does not support <u> underline)
#   - bold horizontal axis text
p <- styled_dotplot(Seurat_object, c("IFIT1B", "VSTM5"), fill_limits = c(0.6, 1.2))
p <- p +
  theme(
    axis.text.x  = element_text(face = "bold"),
    legend.text  = element_markdown(),
    legend.title = element_markdown()
  ) +
  guides(
    fill = guide_legend(
      title        = "<b>Average Expression</b>",
      override.aes = list(shape = 21, colour = "black", stroke = 0.5)
    ),
    size = guide_legend(
      title        = "<b>Percent Expressed</b>",
      override.aes = list(shape = 21, fill = "grey60", colour = "black", stroke = 0.5)
    )
  )
save_pptx(p, "keybiomarkertest.pptx", width = 6.25, height = 5)

# IFIT1/2/3 per cluster 
# Manual steps automated:
#   - shaded percent expressed circles in legend
#   - bold legend category labels
#   - bold horizontal axis text
# NOTE: "Average Expression text nudge up one tick" has no ggplot2 equivalent
#       and remains a minor manual tweak if needed
p <- styled_dotplot(Seurat_object, c("IFIT1", "IFIT2", "IFIT3"), fill_limits = c(2.0, 4.5))
p <- p +
  theme(
    axis.text.x  = element_text(face = "bold"),
    legend.text  = element_markdown(),
    legend.title = element_markdown()
  ) +
  guides(
    fill = guide_legend(
      title        = "<b>Average Expression</b>",
      override.aes = list(shape = 21, colour = "black", stroke = 0.5)
    ),
    size = guide_legend(
      title        = "<b>Percent Expressed</b>",
      override.aes = list(shape = 21, fill = "grey60", colour = "black", stroke = 0.5)
    )
  )
save_pptx(p, "ifittriplegraph.pptx", width = 7, height = 5)

# Top 5 markers per cluster 
# NOTE: manually resize slide in PowerPoint after opening:
#   Design > Slide Size > Custom: width=7, height=7.5 > Ensure Fit
all_markers  <- FindAllMarkers(object = Seurat_object, only.pos = TRUE, min.pct = 0.25)
top5_markers <- all_markers %>% group_by(cluster) %>% slice_head(n = 5) %>% ungroup()
top5_genes   <- unique(top5_markers$gene)
p <- DotPlot(object = Seurat_object, features = top5_genes, scale = FALSE) +
  RotatedAxis() +
  scale_colour_viridis_c(option = "viridis", name = "Average Expression") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = NULL, y = NULL)
doc <- officer::read_pptx()
doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
doc <- officer::ph_with(
  x        = doc,
  value    = rvg::dml(ggobj = p),
  location = officer::ph_location(width = 14, height = 5)
)
print(doc, target = file.path(resultsDir, "longdotplot.pptx"))
message("Saved: ", file.path(resultsDir, "longdotplot.pptx"))


# Cell type per sample 
# NOTE: update sample_name_map with actual sample name mappings for publication
# NOTE: manually shorten sample names and adjust font size in PowerPoint if needed
sample_name_map <- c(
  # "full_sample_name_1" = "Short1"
)
counts_df <- as.data.frame(table(Idents(Seurat_object), Seurat_object@meta.data[["sample"]]))
p <- ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_col(position = position_fill(reverse = FALSE), aes(color = Var1), linewidth = 0.1) +
  labs(x = "Sample", y = "Cell Type Proportion", fill = "Cluster") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    axis.title.x = element_text(margin = margin(t = 15), size = 15),
    axis.title.y = element_text(margin = margin(r = 15), size = 15)
  ) +
  scale_x_discrete(labels = if (length(sample_name_map) > 0) sample_name_map else waiver()) +
  scale_fill_manual(values  = cell_colors, drop = FALSE) +
  scale_color_manual(values = cell_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(color = "none")
save_pptx(p, "bysample.pptx", width = 6.25, height = 4.5)


# Cell type per biomarker category 
# Uses new_ifit (= IFIT1B_Orig_ident) and VSTM5_Orig_ident columns
counts_df <- as.data.frame(table(
  Idents(Seurat_object),
  VSTM5  = Seurat_object@meta.data[["VSTM5_Orig_ident"]],
  IFIT1B = Seurat_object@meta.data[["new_ifit"]]
)) %>%
  pivot_longer(cols = c("VSTM5", "IFIT1B"), names_to = "Category", values_to = "Value") %>%
  mutate(Combination = paste(Category, Value))
counts_summarized <- counts_df %>%
  group_by(Combination, Var1) %>%
  summarize(TotalFreq = sum(Freq), .groups = "drop")
counts_summarized$Var1 <- factor(counts_summarized$Var1, levels = names(cell_colors))
p <- ggplot(counts_summarized, aes(x = Combination, y = TotalFreq, fill = Var1)) +
  geom_col(position = position_fill(reverse = FALSE), aes(color = Var1),
           linewidth = 0.2, linejoin = "mitre") +
  labs(x = NULL, y = "Cell Type Proportion", fill = "Cluster") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.title.x = element_text(margin = margin(t = 15), size = 15),
    axis.title.y = element_text(margin = margin(r = 15), size = 15)
  ) +
  scale_fill_manual(values  = cell_colors, drop = FALSE) +
  scale_color_manual(values = cell_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(color = "none")
save_pptx(p, "bycategory.pptx", width = 6.25, height = 4.5)


# UMAP by cell type 
# NOTE: B cell duplicate label across two sub-clusters must be added manually in PowerPoint
# NOTE: delete legend manually in PowerPoint
p <- DimPlot(Seurat_object, reduction = "umap", cols = cell_colors, label = TRUE) +
  labs(title = "UMAP by Cell Type") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = "none")
save_pptx(p, "umapplot.pptx", width = 8.1, height = 5.9)

# ==========================================================================
# SECTION 2, Plots using CellChat objects
# Replace VSTM5 with IFIT1B throughout to generate equivalent IFIT1B figures
# ==========================================================================

# Outgoing signalling heatmap VSTM5 
# Manual steps automated:
#   - font sizes passed directly to ComplexHeatmap via font.size arguments
#   - title shortened to category name only
# NOTE: drag imported graph fully onto slide in PowerPoint
#       delete top-left "0" on upper scale for spacing
pathway.union <- union(vstm5high@netP$pathways, vstm5low@netP$pathways)
hm1 <- netAnalysis_signalingRole_heatmap(
  vstm5high, pattern = "outgoing", signaling = pathway.union,
  title = "VSTM5 High", font.size = 11, font.size.title = 16
)
hm2 <- netAnalysis_signalingRole_heatmap(
  vstm5low, pattern = "outgoing", signaling = pathway.union,
  title = "VSTM5 Low", font.size = 11, font.size.title = 16
)
save_pptx(hm1 + hm2, "vstm5out.pptx", width = 9.6, height = 5.3, is_complex = TRUE)

# Incoming signalling heatmap VSTM5 
hm3 <- netAnalysis_signalingRole_heatmap(
  vstm5high, pattern = "incoming", signaling = pathway.union,
  title = "VSTM5 High", color.heatmap = "Blues",
  font.size = 11, font.size.title = 16
)
hm4 <- netAnalysis_signalingRole_heatmap(
  vstm5low, pattern = "incoming", signaling = pathway.union,
  title = "VSTM5 Low", color.heatmap = "Blues",
  font.size = 11, font.size.title = 16
)
save_pptx(hm3 + hm4, "vstm5in.pptx", width = 9.6, height = 5.3, is_complex = TRUE)

# CD8 T cell signalling bubble plot VSTM5 
# Manual steps automated:
#   - grey shading placed BEHIND dots by rebuilding ggplot from gg$data
#     (original required "send backwards" manually in PowerPoint)
#   - x axis text size 10, y axis text size 9 italic
#   - p-value legend removed via guides()
#   - "CD8 T Cells ->" label redundancy removed via axis label override
#   - "Commun. Prob." legend title size 12
# NOTE: horizontal axis labels, right-click textbox > Format Shape >
#       Text Options > Text Box > Vertical alignment: Top, then Ctrl+R to align
object.list     <- list(Low = vstm5low, High = vstm5high)
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))
all_clusters    <- levels(cellchat.merged@idents$joint)
cd8_name        <- all_clusters[grep("CD8", all_clusters, ignore.case = TRUE)]

graph_levels <- c(
  "EPCAM\u207a Epithelial", "MKI67\u207a Epithelial", "IFIT\u207a Epithelial",
  "CD8 T Cells", "CD4 T Cells", "Tfh T Cells", "B Cells", "Plasma Cells",
  "Inflammatory monocytes", "Macrophages", "Perivascular Endothelial",
  "Fibroblastic Reticular Cells", "Fibroblasts"
)
cellchat.merged@meta$labels <- factor(
  as.character(cellchat.merged@meta$labels), levels = graph_levels
)
ds_names <- unique(as.character(cellchat.merged@meta$datasets))
cellchat.merged@idents <- list()
cellchat.merged@idents$joint <- cellchat.merged@meta$labels
for (i in ds_names) {
  cellchat.merged@idents[[i]] <- cellchat.merged@meta$labels[
    cellchat.merged@meta$datasets == i
  ]
}

gg <- netVisual_bubble(
  cellchat.merged,
  sources.use = cd8_name,
  targets.use = NULL,
  comparison  = c(2, 1),
  angle.x     = 45,
  thresh      = 0.05
)

if (is.null(gg) || nrow(gg$data) == 0) {
  message("No significant CD8 interactions found, skipping VSTM5CD8out.pptx")
} else {
  message("gg$data columns: ", paste(names(gg$data), collapse=", "))
  
  num_targets <- length(unique(gg$data$target))
  max_limit   <- (num_targets * 2) + 0.5

  # Add shading directly on top of gg (original approach)
  # NOTE: shading sits on top of dots in ggplot, use "send backwards" in PowerPoint
  gg_shaded <- gg +
    annotate(
      "rect",
      xmin  = seq(0.5, max_limit, by = 4),
      xmax  = seq(2.5, max_limit + 2, by = 4),
      ymin  = -Inf, ymax = Inf,
      fill  = "grey92", alpha = 0.5
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.5, max_limit), clip = "on") +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y  = element_text(size = 9, face = "italic"),
      legend.title = element_text(size = 12, face = "bold")
    ) +
    labs(x = NULL, y = NULL)

  save_pptx(gg_shaded, "VSTM5CD8out.pptx", width = 9, height = 5, left = 0.5, top = 1)
}