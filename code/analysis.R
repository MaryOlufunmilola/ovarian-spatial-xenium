# ----------------------------------------------------------
# Xenium Spatial Transcriptomics: Main Analysis
# Authors: Funmi Oyebamiji
# Depends on: load_libraries.R, functions.R 
# Output: /results/Seurat_object.rds
# ----------------------------------------------------------

set.seed(1234)

# Source dependencies if running interactively outside run.sh 
source("/code/load_libraries.R")
source("/code/functions.R")

# ==========================================================================
# STEP 1: Load and pre-process Xenium samples
# ==========================================================================

use_toy_data <- TRUE  # TRUE = Code Ocean toy RDS files
                      # FALSE = real Xenium folders

if (use_toy_data) {
  # Example data path 
  # Expects files named toy_401.rds, toy_402.rds, etc. in /data/
  toy_files  <- list.files(dataDir, pattern = "^toy_\\d+\\.rds$", full.names = TRUE)
  xenium.obj <- list()
  sample_ids <- c()

  for (f in toy_files) {
    j <- gsub("toy_|\\.rds", "", basename(f))   
    message("  Loading toy sample: ", j)
    sample_ids     <- c(sample_ids, j)
    xenium.obj[[j]] <- readRDS(f)
    xenium.obj[[j]]$sample <- j

    # Recalculate nCount and nFeature if missing (toy RDS may not have them)
    if (!"nCount_Xenium" %in% names(xenium.obj[[j]]@meta.data)) {
      xenium.obj[[j]] <- AddMetaData(
        xenium.obj[[j]],
        metadata = Matrix::colSums(xenium.obj[[j]][["Xenium"]]$counts),
        col.name = "nCount_Xenium"
      )
      xenium.obj[[j]] <- AddMetaData(
        xenium.obj[[j]],
        metadata = Matrix::colSums(xenium.obj[[j]][["Xenium"]]$counts > 0),
        col.name = "nFeature_Xenium"
      )
    }

    png(file.path(resultsDir, paste0("VlnPlotbeforefilter_", j, ".png")),
        height = 9, width = 13, units = "in", res = 300)
    print(VlnPlot(xenium.obj[[j]], features = c("nFeature_Xenium", "nCount_Xenium"),
                  layer = "counts", ncol = 2, pt.size = 0, raster = FALSE))
    dev.off()

    print(xenium.obj[[j]])
    saveRDS(xenium.obj[[j]], file.path(resultsDir, paste0(j, ".rds")))
    xenium.obj[[j]] <- subset(xenium.obj[[j]], subset = nCount_Xenium > 20)
  }

} else {
  # Real data path (Xenium folders) 
  xfiles     <- list.files(dataDir)
  xenium.obj <- list()
  sample_ids <- c()

  for (i in xfiles) {
    message("  Preparing spatial data for sample: ", i)

    j <- i 
    sample_ids <- c(sample_ids, j)

    xenium.obj[[j]] <- LoadXenium(paste0(dataDir, i), fov = j, segmentations = "cell")
    xenium.obj[[j]]$sample <- j

    png(file.path(resultsDir, paste0("VlnPlotbeforefilter_", j, ".png")),
        height = 9, width = 13, units = "in", res = 300)
    print(VlnPlot(xenium.obj[[j]], features = c("nFeature_Xenium", "nCount_Xenium"),
                  layer = "counts", ncol = 2, pt.size = 0, raster = FALSE))
    dev.off()

    print(xenium.obj[[j]])
    saveRDS(xenium.obj[[j]], file.path(resultsDir, paste0(j, ".rds")))
    xenium.obj[[j]] <- subset(xenium.obj[[j]], subset = nCount_Xenium > 20)
  }
}

message("Step 1: Done loading samples...")

# ==========================================================================
# STEP 2: Merge, normalise, integrate
# ==========================================================================

xen_integrated <- merge(xenium.obj[[1]], xenium.obj[-1], add.cell.ids = sample_ids)
message("Step 2.1: Done merging samples...")

# Run unique(xen_integrated$sample) after merge to confirm exact names
vstm5_low_samples  <- c("401", "403", "408", "414", "425")  
vstm5_high_samples <- setdiff(unique(xen_integrated$sample), vstm5_low_samples)
ifit1b_low_samples  <- c("401", "402", "403")                          
ifit1b_high_samples <- setdiff(unique(xen_integrated$sample), ifit1b_low_samples)

xen_integrated$VSTM5_Orig_ident  <- ifelse(xen_integrated$sample %in% vstm5_low_samples,  "Low", "High")
xen_integrated$IFIT1B_Orig_ident <- ifelse(xen_integrated$sample %in% ifit1b_low_samples, "Low", "High")

# Normalisation
xen_integrated <- NormalizeData(xen_integrated)
xen_integrated <- FindVariableFeatures(xen_integrated)
xen_integrated <- ScaleData(xen_integrated)
xen_integrated <- SCTransform(xen_integrated,assay = "Xenium")
xen_integrated[["Xenium"]] <- JoinLayers(xen_integrated[["Xenium"]])

# PCA, use variable features selected by SCTransform
xen_integrated <- RunPCA(xen_integrated, npcs = 30,
                          features = VariableFeatures(xen_integrated))
xen_integrated <- RunUMAP(xen_integrated, dims = 1:20, reduction = "pca",
                           reduction.name = "umap.unintegrated")
message("Step 2.2: Done normalising...")

# Harmony integration
xen_integrated <- RunHarmony(xen_integrated, reduction.use = "pca",
                              group.by.vars = "sample", reduction.save = "runharmony")
xen_integrated <- FindNeighbors(xen_integrated, reduction = "runharmony", dims = 1:20)
xen_integrated <- FindClusters(xen_integrated, resolution = 0.03,
                                cluster.name = "runharmony_clusters")
xen_integrated <- RunUMAP(xen_integrated, reduction = "runharmony", dims = 1:20,
                           reduction.name = "umap")
message("Step 2.3: Done integrating...")

#saveRDS(xen_integrated, file.path(resultsDir, "xen_integrated.rds"))

# ==========================================================================
# STEP 3: SingleR cell-type annotation (first-pass, informs manual labels)
# ==========================================================================

# References were pre-downloaded and uploaded to /data/
# To generate locally: saveRDS(celldex::HumanPrimaryCellAtlasData(), "hpca_ref.rds")
hpca.ref   <- readRDS(file.path(dataDir, "hpca_ref.rds"))
blue.ref   <- readRDS(file.path(dataDir, "blueE_ref.rds"))
human_refs <- list("hpca" = hpca.ref, "blueE" = blue.ref)

# Use serial processing for toy data (low cell count), parallel for real data
if (use_toy_data) {
  bpparam <- BiocParallel::SerialParam()
} else {
  bpparam <- BiocParallel::SnowParam(workers = parallelly::availableCores(), type = "SOCK")
}

xen.diet <- as.SingleCellExperiment(xen_integrated, assay = "SCT")

for (ref_name in names(human_refs)) {
  message("Running SingleR with reference: ", ref_name)
  current_ref <- human_refs[[ref_name]]
  cell_main   <- SingleR::SingleR(test = xen.diet, ref = current_ref,
                                   labels = current_ref$label.main, BPPARAM = bpparam)
  xen_integrated[[paste0(ref_name, "_cellmain")]] <- cell_main$labels
}
message("Step 3: Done SingleR annotation...")

# ==========================================================================
# STEP 4: Cluster markers (pre-subclustering)
# ==========================================================================

xen_integrated <- prep_FindMarkers(xen_integrated)
all_markers    <- run_presto(xen_integrated,
                             group_by  = "runharmony_clusters",
                             save_path = file.path(resultsDir, "presto_EachCluster_original.csv"))
xen_integrated$initial_clusters <- xen_integrated$runharmony_clusters
message("Step 4: Done finding cluster markers...")

# ==========================================================================
# STEP 5: Initial Cell type naming (toy data only)
# Real data uses numeric cluster IDs directly in subclustering
# ==========================================================================

if (use_toy_data) {
  # Marker definitions
  canonical_markers <- list(
    "Epithelial" = c("EPCAM", "PPDPF", "SOX17", "CCND1"),
    "Fibroblasts1" = c("DCN", "IGFBP7", "SPARC", "SPARCL1", "ACTA2"),
    "T Cells" = c("CD2", "CD3E", "CD3D", "IL7R", "TIGIT"),
    "Macrophages" = c("CD163", "CD74", "APOE", "FCGR3A", "IFI30"),
    "Plasma Cells" = c("IGHG1", "IGHG3", "IGKC", "MZB1", "JCHAIN"),
    "B Cells" = c("MS4A1", "IGHM", "CD79A", "BANK1", "TNFRSF13C"),
    "Inflammatory monocytes" = c("S100A9", "CXCL8"),
    "Fibroblastic Reticular Cells1" = c("CCL21", "CCL19", "CXCL12")
  )

  # Score modules on downsampled object only
  xen_integrated <- AddModuleScore_UCell(xen_integrated, features = canonical_markers, name = "")

  # Average scores per cluster on subsample
  module_cols <- names(canonical_markers)

  avg_scores <- xen_integrated@meta.data %>%
    group_by(initial_clusters) %>%
    summarise(across(all_of(module_cols), mean), .groups = "drop")

  # Assign label by highest mean score
  cluster_labels_vec <- setNames(
    apply(avg_scores[, module_cols, drop = FALSE], 1,
          function(x) module_cols[which.max(x)]),
    as.character(avg_scores$initial_clusters)
  )

  # Apply to Seurat object 
  xen_integrated@meta.data$Celltypes <- cluster_labels_vec[
    as.character(xen_integrated$initial_clusters)
  ]
  Idents(xen_integrated) <- xen_integrated$Celltypes
  message("Step 5: Done initial cell type naming (toy data)...")
}

# ==========================================================================
# STEP 6: Sub-clustering
# ==========================================================================

if (!use_toy_data) {
  # Real data, subcluster numeric cluster IDs directly
  Idents(xen_integrated) <- xen_integrated$runharmony_clusters

  xen_integrated <- FindSubCluster(xen_integrated, cluster = "0",
                                    graph.name = "SCT_snn",
                                    subcluster.name = "runharmony_subclusters0_005",
                                    resolution = 0.05)
  xen_integrated <- FindSubCluster(xen_integrated, cluster = "1",
                                    graph.name = "SCT_snn",
                                    subcluster.name = "runharmony_subclusters1_002",
                                    resolution = 0.02)
  xen_integrated <- FindSubCluster(xen_integrated, cluster = "3",
                                    graph.name = "SCT_snn",
                                    subcluster.name = "runharmony_subclusters3_002",
                                    resolution = 0.02)

  combined_cluster_labels <- as.character(xen_integrated$runharmony_clusters)
  combined_cluster_labels[xen_integrated$runharmony_clusters == "0"] <-
    as.character(xen_integrated$runharmony_subclusters0_005[xen_integrated$runharmony_clusters == "0"])
  combined_cluster_labels[xen_integrated$runharmony_clusters == "1"] <-
    as.character(xen_integrated$runharmony_subclusters1_002[xen_integrated$runharmony_clusters == "1"])
  combined_cluster_labels[xen_integrated$runharmony_clusters == "3"] <-
    as.character(xen_integrated$runharmony_subclusters3_002[xen_integrated$runharmony_clusters == "3"])

} else {
  # Toy data, subcluster named clusters from Step 5 UCell annotation
  available_clusters <- unique(xen_integrated$Celltypes)
  message("Available clusters: ", paste(available_clusters, collapse=", "))

  if ("Epithelial" %in% available_clusters) {
    xen_integrated <- FindSubCluster(xen_integrated, cluster = "Epithelial",
                                      graph.name = "SCT_snn",
                                      subcluster.name = "runharmony_subclusters0_005",
                                      resolution = 0.05)
  }
  if ("Fibroblasts1" %in% available_clusters) {
    xen_integrated <- FindSubCluster(xen_integrated, cluster = "Fibroblasts1",
                                      graph.name = "SCT_snn",
                                      subcluster.name = "runharmony_subclusters1_002",
                                      resolution = 0.02)
  }
  if ("T Cells" %in% available_clusters) {
    xen_integrated <- FindSubCluster(xen_integrated, cluster = "T Cells",
                                      graph.name = "SCT_snn",
                                      subcluster.name = "runharmony_subclusters3_002",
                                      resolution = 0.02)
  }

  combined_cluster_labels <- as.character(xen_integrated$Celltypes)

  if ("Epithelial" %in% available_clusters) {
    combined_cluster_labels[xen_integrated$Celltypes == "Epithelial"] <-
      as.character(xen_integrated$runharmony_subclusters0_005[xen_integrated$Celltypes == "Epithelial"])
  }
  if ("Fibroblasts1" %in% available_clusters) {
    combined_cluster_labels[xen_integrated$Celltypes == "Fibroblasts1"] <-
      as.character(xen_integrated$runharmony_subclusters1_002[xen_integrated$Celltypes == "Fibroblasts1"])
  }
  if ("T Cells" %in% available_clusters) {
    combined_cluster_labels[xen_integrated$Celltypes == "T Cells"] <-
      as.character(xen_integrated$runharmony_subclusters3_002[xen_integrated$Celltypes == "T Cells"])
  }
}

xen_integrated$runharmony_combined <- combined_cluster_labels

all_markers <- run_presto(xen_integrated,
                          group_by  = "runharmony_combined",
                          save_path = file.path(resultsDir, "presto_EachCluster_final.csv"))
message("Step 5: Done subclustering...")

# ==========================================================================
# STEP 7: Final Cell type naming
# ==========================================================================

Idents(xen_integrated) <- xen_integrated$runharmony_combined

if (!use_toy_data) {
  # Real data, direct renaming matching publication exactly
  new_names <- c(
    "0_0" = "EPCAM\u207A Epithelial", "0_1" = "MKI67\u207A Epithelial",
    "0_2" = "IFIT\u207A Epithelial",  "0_3" = "EPCAM\u207A Epithelial",
    "0_4" = "EPCAM\u207A Epithelial", "0_5" = "EPCAM\u207A Epithelial",
    "1_0" = "Fibroblasts",            "1_1" = "Perivascular Endothelial",
    "1_2" = "Fibroblasts",            "2"   = "Macrophages",
    "3_0" = "CD4 T Cells",            "3_1" = "CD8 T Cells",
    "3_2" = "Tfh T Cells",            "4"   = "Plasma Cells",
    "5"   = "B Cells",                "6"   = "Fibroblastic Reticular Cells",
    "7"   = "Inflammatory Monocytes",
    "8"   = "EPCAM\u207A Epithelial"
  )
  valid_names    <- new_names[names(new_names) %in% levels(Idents(xen_integrated))]
  xen_integrated <- RenameIdents(xen_integrated, valid_names)
  xen_integrated$Celltypes <- as.character(Idents(xen_integrated))
  message("Real data: ", length(valid_names), " clusters renamed")

} else {
  # Toy data, UCell scoring within parent groups
  xen_integrated$InitialCelltypes <- xen_integrated$Celltypes

  subcluster_markers <- list(
    # Epithelial subtypes
    "EPCAM\u207a Epithelial"  = c("EPCAM", "PPDPF", "CCND1"),
    "MKI67\u207a Epithelial"  = c("MKI67", "CDK1", "UBE2C"),
    "IFIT\u207a Epithelial"   = c("IDO1", "CXCL10", "IFIT3", "IRF1"),
    # Fibroblast subtypes
    "Fibroblasts"                  = c("DCN", "C1S", "IGFBP7"),
    "Perivascular Endothelial"     = c("RGS5", "SPARCL1", "FLT1"),
    "Fibroblastic Reticular Cells" = c("CCL21", "CCL19", "CXCL12"),
    # T cell subtypes
    "CD8 T Cells"  = c("CD8A", "GZMA", "GNLY"),
    "CD4 T Cells"  = c("CD4", "IL7R", "LTB"),
    "Tfh T Cells"  = c("CXCL13", "MS4A1", "CORO1A")
  )

  subcluster_module_cols <- names(subcluster_markers)

  xen_integrated <- AddModuleScore_UCell(
    xen_integrated,
    features = subcluster_markers,
    name     = ""
  )

  # average per subcluster
  avg_sub <- xen_integrated@meta.data %>%
    group_by(runharmony_combined) %>%   # replace with your actual subcluster column
    summarise(across(all_of(subcluster_module_cols), mean), .groups = "drop")

  # assign by highest score within parent cell type
  # split by parent so Epithelial subtypes only compete with each other
  epithelial_cols <- c("EPCAM\u207a Epithelial", "MKI67\u207a Epithelial", "IFIT\u207a Epithelial")
  fibroblast_cols <- c("Fibroblasts", "Perivascular Endothelial", "Fibroblastic Reticular Cells")
  tcell_cols      <- c("CD8 T Cells", "CD4 T Cells", "Tfh T Cells")

  assign_within_parent <- function(avg_scores, parent_prefix, cols) {
    avg_scores %>%
      dplyr::filter(grepl(paste0("^", parent_prefix), as.character(runharmony_combined))) %>%
      rowwise() %>%
      mutate(label = cols[which.max(c_across(all_of(cols)))]) %>%
      ungroup() %>%
      dplyr::select(runharmony_combined, label)
  }

  epi_labels   <- if (any(grepl("^Epithelial",  avg_sub$runharmony_combined)))
                    assign_within_parent(avg_sub, "Epithelial",  epithelial_cols) else NULL
  fib_labels   <- if (any(grepl("^Fibroblasts", avg_sub$runharmony_combined)))
                    assign_within_parent(avg_sub, "Fibroblasts", fibroblast_cols) else NULL
  tcell_labels <- if (any(grepl("^T Cells",     avg_sub$runharmony_combined)))
                    assign_within_parent(avg_sub, "T Cells",     tcell_cols)      else NULL

  label_list <- Filter(Negate(is.null), list(epi_labels, fib_labels, tcell_labels))

  if (length(label_list) > 0) {
    all_labels            <- bind_rows(label_list)
    subcluster_labels_vec <- setNames(all_labels$label, all_labels$runharmony_combined)
    message("Subcluster label assignments:")
    print(subcluster_labels_vec)

    xen_integrated@meta.data$Celltypes <- ifelse(
      xen_integrated@meta.data$runharmony_combined %in% names(subcluster_labels_vec),
      subcluster_labels_vec[as.character(xen_integrated@meta.data$runharmony_combined)],
      as.character(xen_integrated@meta.data$runharmony_combined)
    )
  } else {
    message("No subclusters assigned — keeping runharmony_combined as Celltypes")
    xen_integrated@meta.data$Celltypes <- as.character(xen_integrated@meta.data$runharmony_combined)
  }

  # Final NA fallback — runharmony_combined always has values
  na_cells <- is.na(xen_integrated@meta.data$Celltypes)
  if (any(na_cells)) {
    message("Fixing ", sum(na_cells), " NA cells using runharmony_combined")
    xen_integrated@meta.data$Celltypes[na_cells] <-
      as.character(xen_integrated@meta.data$runharmony_combined[na_cells])
  }

  Idents(xen_integrated) <- xen_integrated$Celltypes
}

# Add Celltypes_mod (+ instead of superscript ⁺) for downstream compatibility
#xen_integrated$Celltypes_mod <- gsub("\u207A", "+", xen_integrated$Celltypes)

message("Final Celltypes distribution:")
print(table(xen_integrated$Celltypes, useNA = "ifany"))
message("Step 7: Done cell type naming...")

# ==========================================================================
# STEP 8: Cell type markers
# ==========================================================================

all_markers <- run_presto(xen_integrated,
                          group_by  = "Celltypes",
                          save_path = file.path(resultsDir, "presto_EachCelltype.csv"))

pairwise_markers <- run_presto(xen_integrated,
                               group_by  = "Celltypes",
                               pairwise  = TRUE,
                               save_path = file.path(resultsDir, "presto_EachCelltypevsAnotherCelltype.csv"))

filtered_genes <- all_markers %>%
  dplyr::filter(logFC > 0.25, padj < 0.05)

pathway_gene_list <- filtered_genes$feature %>% unique()

top_genes <- filtered_genes %>%
  dplyr::mutate(group = factor(group, levels = levels(xen_integrated$Celltypes))) %>%
  dplyr::arrange(group, dplyr::desc(logFC)) %>%
  dplyr::group_by(group) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::pull(feature) %>%
  unique()

# Save final annotated object to be loaded by figure_script_1.R and figure_script_2.R
saveRDS(xen_integrated, file.path(resultsDir, "Seurat_object.rds"))
message("Step 7: Done. Seurat_object.rds saved to /results/")
