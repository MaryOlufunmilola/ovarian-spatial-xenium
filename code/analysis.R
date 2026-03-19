# ----------------------------------------------------------
# Xenium Spatial Transcriptomics — Main Analysis
# Author: Funmi Oyebamiji
# Depends on: load_libraries.R, functions.R (sourced by run.sh)
# Output: /results/Seurat_object.rds
# ----------------------------------------------------------

set.seed(1234)

# ── Source dependencies if running interactively outside run.sh ──────────────
if (!exists("dataDir"))    source("/code/load_libraries.R")
if (!exists("save_pptx"))  source("/code/functions.R")

# ==========================================================================
# STEP 1 — Load and pre-process Xenium samples
# ==========================================================================

xfiles     <- list.files(dataDir, pattern = "XETG00270")
xenium.obj <- list()
sample_ids <- c()

for (i in xfiles) {
  message("  Preparing spatial data for sample: ", i)

  j <- gsub("-", ".", strsplit(gsub("INST_1419-", "", i), "__")[[1]][3])
  stopifnot("Sample ID parsing failed — check filename format" = !is.na(j) && nchar(j) > 0)
  sample_ids <- c(sample_ids, j)

  xenium.obj[[j]] <- LoadXenium(paste0(dataDir, i), fov = j, segmentations = "cell")
  xenium.obj[[j]]$sample <- j

  #saveRDS(xenium.obj[[j]], paste0(resultsDir, j, ".rds"))

  # QC filter: remove cells with fewer than 20 transcript counts
  xenium.obj[[j]] <- subset(xenium.obj[[j]], subset = nCount_Xenium > 20)
}
message("Step 1: Done loading samples...")


# ==========================================================================
# STEP 2 — Merge, normalise, integrate
# ==========================================================================

xen_integrated <- merge(xenium.obj[[1]], xenium.obj[-1], add.cell.ids = sample_ids)
message("Step 2.1: Done merging samples...")

# ── Group assignments — replace placeholder names with actual sample j values
# Run unique(xen_integrated$sample) after merge to confirm exact names
vstm5_low_samples  <- c("NM004.00001", "NM004.00003", "NM004.00008", "NM004.00014", "NM004.00025")  
vstm5_high_samples <- setdiff(unique(xen_integrated$sample), vstm5_low_samples)
ifit1b_low_samples  <- c("NM004.00001", "NM004.00002", "NM004.00003")                          
ifit1b_high_samples <- setdiff(unique(xen_integrated$sample), ifit1b_low_samples)

xen_integrated$VSTM5_Orig_ident  <- ifelse(xen_integrated$sample %in% vstm5_low_samples,  "Low", "High")
xen_integrated$IFIT1B_Orig_ident <- ifelse(xen_integrated$sample %in% ifit1b_low_samples, "Low", "High")

# ── Normalisation
xen_integrated <- NormalizeData(xen_integrated)
xen_integrated <- FindVariableFeatures(xen_integrated)
xen_integrated <- ScaleData(xen_integrated)
xen_integrated <- SCTransform(xen_integrated,assay = "Xenium")
xen_integrated[["Xenium"]] <- JoinLayers(xen_integrated[["Xenium"]])

# ── PCA — use variable features selected by SCTransform
xen_integrated <- RunPCA(xen_integrated, npcs = 30,
                          features = VariableFeatures(xen_integrated))
xen_integrated <- RunUMAP(xen_integrated, dims = 1:20, reduction = "pca",
                           reduction.name = "umap.unintegrated")
message("Step 2.2: Done normalising...")

# ── Harmony integration
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
# STEP 3 — SingleR cell-type annotation (first-pass, informs manual labels)
# ==========================================================================

# References must be pre-downloaded and uploaded to /data/
# To generate locally: saveRDS(celldex::HumanPrimaryCellAtlasData(), "hpca_ref.rds")
hpca.ref   <- readRDS(file.path(dataDir, "hpca_ref.rds"))
blue.ref   <- readRDS(file.path(dataDir, "blueE_ref.rds"))
human_refs <- list("hpca" = hpca.ref, "blueE" = blue.ref)

bpparam  <- BiocParallel::SnowParam(workers = parallelly::availableCores(), type = "SOCK")
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
# STEP 4 — Cluster markers (pre-subclustering)
# ==========================================================================

xen_integrated <- prep_FindMarkers(xen_integrated)
all_markers    <- run_presto(xen_integrated,
                             group_by  = "runharmony_clusters",
                             save_path = file.path(resultsDir, "presto_EachCluster_original.csv"))
message("Step 4: Done finding cluster markers...")


# ==========================================================================
# STEP 5 — Sub-clustering
# ==========================================================================

xen_integrated <- FindSubCluster(xen_integrated, cluster = "0", graph.name = "SCT_snn",
                                  subcluster.name = "runharmony_subclusters0_005", resolution = 0.05)
xen_integrated <- FindSubCluster(xen_integrated, cluster = "1", graph.name = "SCT_snn",
                                  subcluster.name = "runharmony_subclusters1_002", resolution = 0.02)
xen_integrated <- FindSubCluster(xen_integrated, cluster = "3", graph.name = "SCT_snn",
                                  subcluster.name = "runharmony_subclusters3_002", resolution = 0.02)

combined_cluster_labels <- as.character(xen_integrated$runharmony_clusters)
combined_cluster_labels[xen_integrated$runharmony_clusters == "0"] <-
  as.character(xen_integrated$runharmony_subclusters0_005[xen_integrated$runharmony_clusters == "0"])
combined_cluster_labels[xen_integrated$runharmony_clusters == "1"] <-
  as.character(xen_integrated$runharmony_subclusters1_002[xen_integrated$runharmony_clusters == "1"])
combined_cluster_labels[xen_integrated$runharmony_clusters == "3"] <-
  as.character(xen_integrated$runharmony_subclusters3_002[xen_integrated$runharmony_clusters == "3"])
xen_integrated$runharmony_combined <- combined_cluster_labels

all_markers <- run_presto(xen_integrated,
                          group_by  = "runharmony_combined",
                          save_path = file.path(resultsDir, "presto_EachCluster_final.csv"))
message("Step 5: Done subclustering...")


# ==========================================================================
# STEP 6 — Cell type naming
# ==========================================================================

Idents(xen_integrated) <- xen_integrated$runharmony_combined

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
xen_integrated$Celltypes <- Idents(xen_integrated)
#xen_integrated$cell_ids <- colnames(xen_integrated)
xen_integrated$Celltypes_mod <- gsub("\u207A", "+", xen_integrated$Celltypes)
message("Step 6: Done cell type naming...")


# ==========================================================================
# STEP 7 — Cell type markers
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

# Save final annotated object — loaded by figure_script_1.R and figure_script_2.R
saveRDS(xen_integrated, file.path(resultsDir, "Seurat_object.rds"))
message("Step 7: Done. Seurat_object.rds saved to /results/")
