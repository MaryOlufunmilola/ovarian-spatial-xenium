# ----------------------------------------------------------
# Shared Utility Functions
# Authors: Funmi Oyebamiji
# Sourced by: analysis.R, figure_script_1.R
# ----------------------------------------------------------
# Usage in each downstream script:
#   source("/code/functions.R")
# ----------------------------------------------------------

# ==========================================================================
# SECTION 1 — Colour palettes
# ==========================================================================

sample_colors <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
  "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
)

cluster_colors <- c(
  "0"   = "#1F78B4", "0_0" = "#084594", "0_1" = "#2171B5",
  "0_2" = "#4292C6", "0_3" = "#6BAED6", "0_4" = "#9ECAE1",
  "0_5" = "#8DA0CB", "1"   = "#06613f", "1_0" = "#33A02C",
  "1_1" = "#A6D854", "1_2" = "#06613f", "2"   = "#E31A1C",
  "3"   = "#FF7F00", "3_0" = "#7F2704", "3_1" = "#A63603",
  "3_2" = "#D94801", "3_3" = "#F16913", "3_4" = "#FD8D3C",
  "4"   = "#6A3D9A", "5"   = "#B15928", "6"   = "#66C2A5",
  "7"   = "#E78AC3", "8"   = "#FB9A99"
)

cell_colors <- c(
  "EPCAM\u207A Epithelial"          = "#CFBEB4",  # ⁺ unicode
  "MKI67\u207A Epithelial"          = "#B1B6BA",
  "IFIT\u207A Epithelial"           = "#968F9F",
  "EPCAM+ Epithelial"               = "#CFBEB4",  # plain + fallback
  "MKI67+ Epithelial"               = "#B1B6BA",
  "IFIT+ Epithelial"                = "#968F9F",
  "Fibroblasts"                     = "#33A02C",
  "Perivascular Endothelial"        = "#A6D854",
  "Macrophages"                     = "#E31A1C",
  "CD4 T Cells"                     = "#A63603",
  "CD8 T Cells"                     = "#FDBF6F",
  "Tfh T Cells"                     = "#FFD92F",
  "B Cells"                         = "#4292C6",
  "Plasma Cells"                    = "#9ECAE1",
  "Fibroblastic Reticular Cells"    = "#66C2A5",
  "Inflammatory Monocytes"          = "#FB9A99"
)

singler_colors <- c(
  RColorBrewer::brewer.pal(9,  "Set1"),
  RColorBrewer::brewer.pal(8,  "Set2"),
  RColorBrewer::brewer.pal(8,  "Dark2"),
  RColorBrewer::brewer.pal(12, "Set3"),
  RColorBrewer::brewer.pal(9,  "Pastel1")
)


# ==========================================================================
# SECTION 2 — Differential expression
# ==========================================================================

#' Prepare SCTransform model slots for FindMarkers / run_presto
#'
#' Sets the umi.assay slot on each SCTModel to "Xenium" then runs
#' PrepSCTFindMarkers(). Required before any marker-finding step when
#' SCTransform was run on a Xenium assay.
#'
#' @param obj  Seurat object with SCT assay
#' @return     Seurat object ready for FindMarkers / wilcoxauc
prep_FindMarkers <- function(obj) {
  num_slices <- length(obj@assays$SCT@SCTModel.list)
  for (i in 2:num_slices) {
    obj@assays$SCT@SCTModel.list[[i]]@umi.assay <- "Xenium"
  }
  obj <- obj %>% PrepSCTFindMarkers()
  return(obj)
}

#' Run Presto Wilcoxon test with optional pairwise comparisons
#'
#' @param seurat_obj     Seurat object
#' @param group_by       Metadata column to group by
#' @param assay          Assay to use (default "SCT")
#' @param groups_use     Subset of groups to test (default NULL = all)
#' @param add_pct_diff   Add percent difference column (default TRUE)
#' @param save_path      Full path to save CSV output — use /results/filename.csv
#' @param pairwise       Run all pairwise comparisons (default FALSE)
run_presto <- function(seurat_obj,
                       group_by,
                       assay        = "SCT",
                       groups_use   = NULL,
                       add_pct_diff = TRUE,
                       save_path    = NULL,
                       pairwise     = FALSE) {
  if (pairwise) {
    clusters      <- unique(seurat_obj[[group_by]][, 1])
    cluster_pairs <- combn(clusters, 2, simplify = FALSE)
    results <- purrr::map_dfr(cluster_pairs, function(pair) {
      res <- presto::wilcoxauc(seurat_obj, group_by = group_by,
                               seurat_assay = assay, groups_use = pair)
      if (add_pct_diff) res <- scCustomize::Add_Pct_Diff(res, pct.1_name = "pct_in", pct.2_name = "pct_out")
      res$cluster_1 <- pair[1]
      res$cluster_2 <- pair[2]
      res
    })
  } else {
    results <- presto::wilcoxauc(seurat_obj, group_by = group_by,
                                 seurat_assay = assay, groups_use = groups_use)
    if (add_pct_diff) results <- scCustomize::Add_Pct_Diff(results, pct.1_name = "pct_in", pct.2_name = "pct_out")
  }
  if (!is.null(save_path)) {
    write.csv(results, file = save_path, row.names = FALSE)
    message("Saved: ", save_path)
  }
  return(results)
}


# ==========================================================================
# SECTION 3 — Plot savers
# ==========================================================================

#' Save a DimPlot to PNG
#'
#' @param seurat_obj  Seurat object
#' @param filename    Full output path — use /results/filename.png
#' @param reduction   Reduction to use e.g. "umap"
#' @param group_by    Metadata column to colour by
#' @param ...         Additional DimPlot arguments
#' @param width       PNG width in inches (default 13)
#' @param height      PNG height in inches (default 9)
save_dimplot <- function(seurat_obj, filename, reduction, group_by,
                         split_by = NULL, raster = FALSE, cols,
                         shape_by = NULL, order = NULL, shuffle = FALSE,
                         label = FALSE, na_value = NULL, label_size = 4,
                         ncol = NULL, pt_size = 0.001,
                         cells_highlight = NULL, cols_highlight = NULL,
                         sizes_highlight = NULL, width = 13, height = 9) {
  if (is.null(na_value)) na_value <- "lightgrey"
  plot_obj <- Seurat::DimPlot(
    seurat_obj,
    reduction        = reduction,
    group.by         = group_by,
    split.by         = split_by,
    shape.by         = shape_by,
    label            = label,
    label.size       = label_size,
    ncol             = ncol,
    pt.size          = pt_size,
    order            = order,
    shuffle          = shuffle,
    raster           = raster,
    cells.highlight  = cells_highlight,
    cols.highlight   = cols_highlight,
    sizes.highlight  = sizes_highlight,
    na.value         = na_value,
    cols             = cols
  )
  png(filename, height = height, width = width, units = "in", res = 300)
  print(plot_obj)
  dev.off()
  message("Saved: ", filename)
}

#' Save a scCustomize DotPlot to PNG
#'
#' @param seurat_obj  Seurat object
#' @param filename    Full output path — use /results/filename.png
#' @param features    Gene features to plot
#' @param group_by    Metadata column to group by
#' @param width       PNG width in inches (default 13)
#' @param height      PNG height in inches (default 9)
save_dotplot <- function(seurat_obj, filename, width = 13, height = 9,
                         features, group_by, x_lab_rotate, x_size, y_size,
                         colors = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(100)) {
  png(filename, height = height, width = width, units = "in", res = 300)
  print(
    scCustomize::DotPlot_scCustom(
      seurat_obj,
      features    = features,
      group.by    = group_by,
      colors_use  = colors,
      x_lab_rotate = x_lab_rotate
    ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = x_size),
        axis.text.y = ggplot2::element_text(size = y_size)
      )
  )
  dev.off()
  message("Saved: ", filename)
}

#' Save a DoHeatmap to PNG using a downsampled set of cells
#'
#' @param seurat_obj    Seurat object
#' @param features      Gene features for heatmap rows
#' @param filename      Full output path — use /results/filename.png
#' @param group_by      Metadata column to group cells (default "runharmony_celltypes")
#' @param group_colors  Named colour vector for groups
#' @param numcells      Cells to downsample per group (default 100)
#' @param width         PNG width in inches (default 13)
#' @param height        PNG height in inches (default 9)
#' @param ysize         Y-axis label font size (default 5)
save_heatmap <- function(seurat_obj, features, filename,
                         group_by     = "runharmony_celltypes",
                         group_colors = NULL,
                         numcells     = 100,
                         width        = 13,
                         height       = 9,
                         ysize        = 5) {
  plot_obj <- Seurat::DoHeatmap(
    seurat_obj,
    features     = features,
    group.by     = group_by,
    cells        = scCustomize::Random_Cells_Downsample(seurat_obj, num_cells = numcells, allow_lower = TRUE),
    group.colors = group_colors,
    raster       = FALSE
  ) +
    ggplot2::scale_fill_gradientn(
      colors = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
    ) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = ysize))
  png(filename, height = height, width = width, units = "in", res = 300)
  print(plot_obj)
  dev.off()
  message("Saved: ", filename)
}


# ==========================================================================
# SECTION 4 — CellChat pipeline
# ==========================================================================

#' Run the full CellChat spatial pipeline on a named list of Seurat objects
#'
#' Image-based tissue coordinates merged across all images,
#' sample column renamed to 'samples', cell type labels from 'Celltypes' column.
#'
#' @param seurat_list       Named list of Seurat objects (one per condition)
#' @param celltype_col      Metadata column for cell type labels (default "Celltypes")
#' @param assay             Assay to pull expression data from (default "Xenium").
#'                          Change to "SCT" to use SCTransform counts instead.
#' @param contact.range     CellChat contact range (default 100)
#' @param interaction.range CellChat interaction range (default 250)
#' @param scale.distance    CellChat scale distance (default 7.4)
#' @param nboot             Number of bootstraps for computeCommunProb (default 100)
#' @param workers           Number of cores for parallel processing. Defaults to
#'                          parallelly::availableCores() — automatically detects
#'                          cores available in the current environment (local or
#'                          Code Ocean). Override with e.g. workers = 8 if needed.
#' @return Named list of CellChat objects
run_cellchat_pipeline <- function(seurat_list,
                                  celltype_col      = "Celltypes",
                                  assay             = "Xenium",
                                  contact.range     = 100,
                                  interaction.range = 250,
                                  scale.distance    = 7.4,
                                  nboot             = 100,
                                  workers           = parallelly::availableCores()) {
  cellchat_list <- lapply(names(seurat_list), function(cond) {
    obj <- seurat_list[[cond]]
    message("Processing condition: ", cond)

    # ── Pull data ──────────────────────────────────────────────────────
    data.input <- Seurat::GetAssayData(obj, slot = "data", assay = assay)

    # ── Merge tissue coordinates across all images ─────────────────────────
    image_names <- Seurat::Images(obj)
    all_tissue_coordinates <- setNames(
      lapply(image_names, function(img_name) {
        Seurat::GetTissueCoordinates(object = obj, image = img_name)
      }),
      image_names
    )
    merged_cords           <- data.table::rbindlist(all_tissue_coordinates)
    rownames(merged_cords) <- merged_cords$cell
    merged_cords$cell      <- NULL

    # ── Prepare metadata ───────────────────────────────────────────────────
    meta <- obj@meta.data
    colnames(meta)[colnames(meta) == "sample"] <- "samples"
    meta$samples <- as.factor(meta$samples)
    meta$labels  <- meta[[celltype_col]]

    # ── Create CellChat object ─────────────────────────────────────────────
    spatial.factors <- data.frame(ratio = 1, tol = 10)
    cellchat <- CellChat::createCellChat(
      object          = data.input,
      meta            = meta,
      group.by        = "labels",
      datatype        = "spatial",
      coordinates     = merged_cords,
      spatial.factors = spatial.factors
    )

    # ── Database and gene filtering ────────────────────────────────────────
    cellchat@DB <- CellChat::subsetDB(CellChat::CellChatDB.human)
    cellchat    <- CellChat::subsetData(cellchat, features = NULL)

    future::plan("multisession", workers = workers)
    cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
    cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

    # ── Communication probability ──────────────────────────────────────────
    # Range 250 parameters (default). For range 50 use:
    #   contact.range=15, interaction.range=50, scale.distance=0.5
    cellchat <- CellChat::computeCommunProb(
      cellchat,
      type              = "truncatedMean",
      trim              = 0.1,
      distance.use      = TRUE,
      contact.range     = contact.range,
      interaction.range = interaction.range,
      contact.dependent = TRUE,
      scale.distance    = scale.distance,
      nboot             = nboot
    )
    future::plan("sequential")

    # ── Post-processing ────────────────────────────────────────────────────
    cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
    cellchat <- CellChat::computeCommunProbPathway(cellchat)
    cellchat <- CellChat::aggregateNet(cellchat)
    cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")

    message("Done: ", cond)
    return(cellchat)
  })
  names(cellchat_list) <- names(seurat_list)
  return(cellchat_list)
}

#' Generate and save CellChat heatmaps and bubble plots for a condition pair
#'
#' @param rds_path    Full path to saved CellChat list RDS — use /results/filename.rds
#' @param prefix      Output file prefix — PNGs saved to /results/prefix_*.png
#' @param cell_colors Named colour vector for cell types
analyze_cellchat <- function(rds_path, prefix, cell_colors) {
  cellchat_list <- readRDS(rds_path)
  i             <- 1
  pathway.union <- union(
    cellchat_list[[i]]@netP$pathways,
    cellchat_list[[i + 1]]@netP$pathways
  )

  ht_outgoing <- lapply(1:2, function(j) {
    CellChat::netAnalysis_signalingRole_heatmap(
      cellchat_list[[i + j - 1]],
      pattern    = "outgoing",
      signaling  = pathway.union,
      title      = names(cellchat_list)[i + j - 1],
      width      = 12, height = 9,
      color.use  = cell_colors[rownames(cellchat_list[[i]]@netP$prob)]
    )
  })

  ht_incoming <- lapply(1:2, function(j) {
    CellChat::netAnalysis_signalingRole_heatmap(
      cellchat_list[[i + j - 1]],
      pattern       = "incoming",
      signaling     = pathway.union,
      title         = names(cellchat_list)[i + j - 1],
      width         = 12, height = 9,
      color.heatmap = "GnBu",
      color.use     = cell_colors[rownames(cellchat_list[[i]]@netP$prob)]
    )
  })

  # Save heatmaps to /results/
  png(file.path("/results", paste0(prefix, "_outgoing_heatmaps.png")),
      width = 13, height = 9, res = 300, units = "in")
  ComplexHeatmap::draw(ht_outgoing[[1]] + ht_outgoing[[2]],
                       ht_gap = grid::unit(0.75, "cm"))
  dev.off()

  png(file.path("/results", paste0(prefix, "_incoming_heatmaps.png")),
      width = 13, height = 9, res = 300, units = "in")
  ComplexHeatmap::draw(ht_incoming[[1]] + ht_incoming[[2]],
                       ht_gap = grid::unit(0.75, "cm"))
  dev.off()

  # Bubble plots
  cellchat <- CellChat::mergeCellChat(cellchat_list, add.names = names(cellchat_list))

  gg <- CellChat::netVisual_bubble(
    cellchat,
    comparison  = c(1, 2),
    sources.use = c(3),
    targets.use = c(1:4, 7:10, 12:13),
    angle.x     = 45
  )
  png(file.path("/results", paste0(prefix, "_signaling_heatmaps_general.png")),
      width = 13, height = 9, res = 300, units = "in")
  print(gg)
  dev.off()

  gg1 <- CellChat::netVisual_bubble(
    cellchat,
    comparison   = c(1, 2),
    sources.use  = c(1:3, 5, 7:9, 12, 13),
    targets.use  = c(2, 3, 13),
    max.dataset  = 2,
    title.name   = paste("Increased signaling in", prefix, "High"),
    angle.x      = 45,
    remove.isolate = TRUE
  )
  gg2 <- CellChat::netVisual_bubble(
    cellchat,
    comparison   = c(1, 2),
    sources.use  = c(1:3, 5, 7:9, 12, 13),
    targets.use  = c(2, 3, 13),
    max.dataset  = 1,
    title.name   = paste("Decreased signaling in", prefix, "High"),
    angle.x      = 45,
    remove.isolate = TRUE
  )
  png(file.path("/results", paste0(prefix, "_signaling_heatmaps_diff.png")),
      width = 13, height = 11, res = 300, units = "in")
  print(gg1 / gg2)
  dev.off()

  message("analyze_cellchat complete — outputs saved to /results/", prefix, "_*.png")
}


# ==========================================================================
# SECTION 5 — Cell type proportion analysis
# ==========================================================================

#' Compute, save, and plot cell type proportions
#'
#' @param data          Data frame (typically Seurat object metadata)
#' @param group_vars    Character vector of grouping column names
#' @param fill_var      Column name for fill (cell type)
#' @param plot_title    Title string for the plot
#' @param output_prefix File prefix — CSV and PNG saved to /results/prefix.*
#' @param facet_var     Optional column name for faceting
analyze_proportions <- function(data, group_vars, fill_var, plot_title,
                                output_prefix, facet_var = NULL) {
  prop_df <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars)), !!rlang::sym(fill_var)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::mutate(proportion = count / sum(count)) %>%
    dplyr::ungroup()

  csv_path <- file.path("/results", paste0(output_prefix, ".csv"))
  write.csv(prop_df, csv_path, row.names = FALSE)
  message("Saved: ", csv_path)

  p <- ggplot2::ggplot(
    prop_df,
    ggplot2::aes_string(x = tail(group_vars, 1), y = "proportion", fill = fill_var)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::labs(y = "Cell Type Proportion", x = "", fill = "", title = plot_title) +
    ggplot2::scale_fill_manual(values = cell_colors) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text.x  = ggplot2::element_text(size = 14, angle = 45, hjust = 0.75),
      axis.text.y  = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      panel.grid   = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = 18),
      legend.title = ggplot2::element_blank()
    )

  if (!is.null(facet_var)) {
    p <- p +
      ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = "free_x") +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 14, face = "bold"))
  }

  png_path <- file.path("/results", paste0(output_prefix, ".png"))
  ggplot2::ggsave(png_path, p, width = 13, height = 9, dpi = 300)
  message("Saved: ", png_path)
}


# ==========================================================================
# SECTION 6 — Cluster composition table and plot
# ==========================================================================

#' Build a wide-format cluster composition table from a Seurat object
#'
#' @param seurat_object  Seurat object
#' @param cluster_column Metadata column for clusters/cell types
#' @param sample_column  Metadata column for sample identity
#' @return Wide data frame: rows = samples, columns = clusters + total_cell_count
create_cluster_composition_table <- function(seurat_object,
                                             cluster_column,
                                             sample_column) {
  table_out <- seurat_object@meta.data %>%
    dplyr::group_by(!!rlang::sym(sample_column), !!rlang::sym(cluster_column)) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      id_cols    = !!rlang::sym(sample_column),
      names_from = !!rlang::sym(cluster_column),
      values_from = "Count"
    ) %>%
    dplyr::left_join(
      seurat_object@meta.data %>%
        dplyr::group_by(!!rlang::sym(sample_column)) %>%
        dplyr::summarise(total_cell_count = dplyr::n(), .groups = "drop"),
      by = sample_column
    )

  colnames(table_out)[-c(1, ncol(table_out))] <-
    paste0("", sort(unique(seurat_object@meta.data[[cluster_column]])))

  return(table_out)
}

#' Plot cluster composition as a stacked bar chart
#'
#' @param table_samples_by_clusters Output from create_cluster_composition_table()
#' @param sample_column   Column name for samples (x axis)
#' @param plotTitle       Plot title string
#' @param cluster_column  Label for the fill legend
#' @param cluster_order   Optional factor order for clusters
#' @param type            "proportion" or "count"
#' @param cluster_colors  Named colour vector
#' @param save_plot       Save to file if TRUE
#' @param filename        Full output path — use /results/filename.png
#' @param width           PNG width in inches (default 13)
#' @param height          PNG height in inches (default 9)
#' @param dpi             PNG resolution (default 300)
plot_cluster_composition <- function(table_samples_by_clusters,
                                     sample_column,
                                     plotTitle       = "title",
                                     cluster_column  = "Cluster",
                                     cluster_order   = NULL,
                                     type            = "count",
                                     cluster_colors  = NULL,
                                     save_plot       = FALSE,
                                     filename        = NULL,
                                     width           = 13,
                                     height          = 9,
                                     dpi             = 300) {
  cluster_cols <- setdiff(
    colnames(table_samples_by_clusters),
    c(sample_column, "total_cell_count")
  )
  .data <- table_samples_by_clusters

  if (type == "proportion") {
    .data[, cluster_cols] <- .data[, cluster_cols] / .data$total_cell_count
    value_column <- "Proportion"
  } else {
    value_column <- "Count"
  }

  .data <- tidyr::pivot_longer(
    .data,
    cols      = tidyr::all_of(cluster_cols),
    names_to  = cluster_column,
    values_to = value_column
  )

  if (!is.null(cluster_order)) {
    .data[[cluster_column]] <- factor(.data[[cluster_column]], levels = cluster_order)
  }

  p <- if (type == "proportion") {
    ggplot2::ggplot(
      .data,
      ggplot2::aes(x = .data[[sample_column]], y = .data[[value_column]],
                   fill = .data[[cluster_column]])
    ) +
      ggplot2::geom_bar(stat = "identity", position = "fill") +
      ggplot2::theme_classic() +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(title = plotTitle, x = "", y = "Proportion", fill = cluster_column)
  } else {
    ggplot2::ggplot(
      .data,
      ggplot2::aes(x = .data[[sample_column]], y = .data[[value_column]],
                   fill = .data[[cluster_column]])
    ) +
      ggplot2::geom_bar(stat = "identity", position = "stack") +
      ggplot2::theme_classic() +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(title = plotTitle, x = "", y = "Cell Count", fill = cluster_column)
  }

  if (!is.null(cluster_colors)) {
    p <- p + ggplot2::scale_fill_manual(values = cluster_colors)
  }

  if (save_plot && !is.null(filename)) {
    ggplot2::ggsave(filename, plot = p, width = width, height = height, dpi = dpi)
    message("Saved: ", filename)
  }

  return(p)
}
