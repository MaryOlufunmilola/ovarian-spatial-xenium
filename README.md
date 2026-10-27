# Xenium Spatial Transcriptomics: Olaparib plus Tremelimumab in Recurrent Ovarian Cancer

R code for the Xenium spatial transcriptomics analysis in:

> Adams SF, Gimotty PA, Rutledge TL, et al. Olaparib plus tremelimumab in women with
> recurrent ovarian cancer and a germline *BRCA1* or *BRCA2* mutation: a phase I/II study.
> *Nature Communications* (2026). https://doi.org/10.1038/s41467-026-78395-y

Code by Funmi Oyebamiji and Eleanor Paskus.

The paper reports, as an exploratory post hoc analysis, that low VSTM5 or IFIT1B
expression was associated with treatment benefit and longer progression free survival.
This repository covers the spatial part of that work: loading and QC of 8 Xenium
samples, integration, clustering, cell type annotation, spatial cell to cell
communication (CellChat) in the VSTM5 low and high groups, and the figures.

## What the code does

| Step | Script | What it does |
|---|---|---|
| 1 | `analysis.R` | Loads each sample, plots counts and genes per cell, keeps cells with more than 20 counts |
| 2 | `analysis.R` | Merges samples, assigns each sample to the VSTM5 and IFIT1B low or high group, SCTransform, PCA, Harmony integration across samples, Leiden clustering, UMAP |
| 3 | `analysis.R` | First pass labels from SingleR (Human Primary Cell Atlas and Blueprint/ENCODE references) |
| 4 to 7 | `analysis.R` | Cluster markers (presto), subclustering of epithelial, fibroblast, and T cell clusters, final cell type names |
| 8 | `analysis.R` | Markers for each cell type and for each pair of cell types; saves `Seurat_object.rds` |
| 9 | `figure_script_1.R` | Spatial CellChat, run separately for the VSTM5 low and high groups |
| 10 | `figure_script_2.R` | Paper figures (dot plots, cell type composition, CellChat comparisons), saved as PowerPoint files |
| optional | `eni_figures.R` | Diagnostic UMAPs (by sample, cluster, and SingleR label); not in the paper |

`run` runs steps 1 to 10 in order and records each step's run time in `timing.tsv` in
the results folder.

Cell types in the paper were named by inspecting marker genes for each cluster and
subcluster (`new_names` in Step 7 of `analysis.R`). SingleR labels were used as a guide,
not as the final labels.

Spatial CellChat uses the cell coordinates with `interaction.range = 250` and
`scale.distance = 7.4` (see `run_cellchat_pipeline()` in `functions.R`).

## Running it

### On Code Ocean

The analysis is packaged as a Code Ocean capsule with the environment in
`environment/Dockerfile`. A **Reproducible Run** calls `run`, which uses `/code`,
`/data`, and `/results`.

The capsule runs on a small toy subset of the data (files named `toy_401.rds`,
`toy_402.rds`, and so on) so that it fits in the capsule's memory. On the toy data,
cell types are assigned by UCell scores for marker sets, because the manual cluster
names in Step 7 only apply to the full dataset. Toy results show that the code runs
end to end; they will not match the paper's figures.

### On your own machine

You need the Xenium output folders (see Data below), R 4.4, and the packages in
`environment/Dockerfile`. The full dataset needs more than 128 GB of RAM.

```bash
git clone https://github.com/MaryOlufunmilola/ovarian-spatial-xenium.git
cd ovarian-spatial-xenium
mkdir -p data results
# put one Xenium output folder per sample in data/
```

Save the SingleR references to `data/` once:

```r
saveRDS(celldex::HumanPrimaryCellAtlasData(), "data/hpca_ref.rds")
saveRDS(celldex::BlueprintEncodeData(),       "data/blueE_ref.rds")
```

Then run everything from the repository root:

```bash
CODE_DIR=code DATA_DIR=data/ RESULTS_DIR=results/ USE_TOY_DATA=FALSE bash run
```

Without these settings the scripts use the Code Ocean paths and the toy data, so the
capsule runs unchanged.

## Inputs and outputs

| | Files |
|---|---|
| Input | One Xenium output folder per sample (or `toy_<id>.rds` files); `hpca_ref.rds` and `blueE_ref.rds` |
| Main object | `Seurat_object.rds`: integrated object with clusters, SingleR labels, final `Celltypes`, and the VSTM5 and IFIT1B groups |
| CellChat | `vstm5low.rds`, `vstm5high.rds` |
| Markers | `presto_EachCluster_original.csv`, `presto_EachCluster_final.csv`, `presto_EachCelltype.csv`, `presto_EachCelltypevsAnotherCelltype.csv` |
| Figures | QC violin plots for each sample (PNG); paper figures as PowerPoint files from `figure_script_2.R` |
| Run record | `timing.tsv` |

## Software

The environment is defined in `environment/Dockerfile`, built on Code Ocean's R 4.4.2
image. Most CRAN packages are pinned to a version (for example Seurat 5.5.1,
scCustomize 3.3.0), and harmony, presto, ComplexHeatmap, and circlize are pinned to
GitHub commits. SeuratObject, tidyverse, officer, rvg, NMF, CellChat, and the
Bioconductor packages (SingleR, UCell, glmGamPoi, BiocParallel) install the version
current at build time, so a rebuild can pick up newer versions of those.

## Data

Xenium data from the trial's patient samples are not in this repository. See the Data
availability statement in the paper for how to access them.

## Structure

```
code/load_libraries.R   packages, paths, options
code/functions.R        shared functions (markers, plots, CellChat, proportions)
code/analysis.R         steps 1 to 8
code/figure_script_1.R  spatial CellChat for one VSTM5 group
code/figure_script_2.R  paper figures (written by Eleanor Paskus, modified by Funmi Oyebamiji)
code/eni_figures.R      optional diagnostic figures
environment/Dockerfile  Code Ocean environment
run                     runs everything in order
```

## Citation

If you use this code, please cite the paper above.

## License

MIT; see `LICENSE`.
