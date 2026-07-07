# Spatial Transcriptomics Analysis of Ovarian Cancer Patient Samples

> **Manuscript title:** [PLACEHOLDER — update when accepted]  
> **Journal:** Nature Communications (submitted)  
> **Code authors:** Funmi Oyebamiji, Eleanor Paskus

---

## Overview

This repository contains the R code for the Xenium spatial transcriptomics component of a multi-modal study of ovarian cancer patient samples. The codebase covers spatial data loading, quality control, integration, clustering, cell-type annotation, cell-cell communication analysis (CellChat), and figure generation.

A fully reproducible computational capsule is available on Code Ocean:  
**Code Ocean DOI:** [add after capsule is published]

---

## Repository structure

```
.
├── code/
│   ├── load_libraries.R      # Package loading and path setup
│   ├── functions.R           # Shared utility functions (run_presto, run_cellchat_pipeline, etc.)
│   ├── analysis.R            # Main analysis: QC, integration, clustering, annotation
│   ├── figure_script_1.R     # CellChat spatial pipeline per condition 
│   └── figure_script_2.R     # Manuscript figures from Seurat and CellChat objects (written by Eleanor Paskus, modified by Funmi Oyebamiji)
├── environment/
│   └── Dockerfile            # Reproducible R environment for Code Ocean
├── run.sh                    # Code Ocean entry point, executes scripts in order
└── README.md
```

> **Note on data:** Raw and processed Xenium files are not hosted here due to file size and patient privacy constraints. See the Data Availability section below.

> **Note on toy data:** The Code Ocean capsule runs on a small representative subset of the data (toy dataset) due to RAM constraints. Full results require running on the complete dataset locally or on a high-memory system (>128 GB RAM).

---

## Requirements

### R version
R ≥ 4.4.2 (tested on R 4.4.2 in the Code Ocean environment)

### Key dependencies

| Package | Source | Version | Purpose |
|---|---|---|---|
| Seurat | CRAN | 5.5.1 | Spatial object handling, clustering, SCTransform |
| SeuratObject | CRAN | latest | Seurat v5 object structure |
| scCustomize | CRAN | 3.3.0 | Seurat visualisation helpers |
| SingleR | Bioconductor | 3.20 | Automated cell-type annotation |
| UCell | Bioconductor | 3.20 | Gene signature scoring |
| harmony | GitHub | pinned commit | Batch integration |
| presto | GitHub | pinned commit | Fast Wilcoxon marker testing |
| CellChat | GitHub (jinworks) | latest | Cell-cell communication |
| ComplexHeatmap | GitHub (jokergoo) | pinned commit | Heatmap figures |
| circlize | GitHub (jokergoo) | pinned commit | Circular visualisation |
| BiocParallel | Bioconductor | 3.20 | Parallel processing |
| glmGamPoi | Bioconductor | 3.20 | Fast SCTransform v2 estimation |
| NMF | CRAN | latest | Required by CellChat |
| officer | CRAN | latest | PowerPoint export |
| rvg | CRAN | latest | Vector graphics for PowerPoint |
| ggtext | CRAN | latest | Markdown in ggplot2 labels |
| tidyverse | CRAN | latest | Data wrangling and plotting |

All packages and their exact versions are managed via the `Dockerfile` in `environment/`. GitHub-sourced packages are pinned to specific commit hashes for full reproducibility.

---

## Reproducing the analysis

### On Code Ocean (recommended)

1. Open the Code Ocean capsule at the DOI above
2. Xenium toy data files are pre-loaded under `/data/`
3. Reference files (`hpca_ref.rds`, `blueE_ref.rds`) for SingleR are pre-loaded under `/data/`
4. Click **Reproducible Run** — scripts execute automatically via `run.sh` in this order:

```bash
Rscript /code/analysis.R
Rscript /code/figure_script_1.R Low
Rscript /code/figure_script_1.R High
Rscript /code/figure_script_2.R
```

> **RAM requirement:** The full pipeline requires >128 GB RAM on real data. On Code Ocean, the toy dataset runs within the available memory. `figure_script_1.R` is run sequentially per condition (Low then High) to avoid simultaneous memory usage.

> **Note on toy data results:** Cell type labels on toy data are assigned via UCell gene signature scoring (an approximation). On the full dataset, cell types are assigned by direct cluster renaming based on validated marker gene inspection. Results on toy data are therefore indicative rather than identical to the published figures.

---

### Locally (full dataset)

1. Clone this repository:

```bash
git clone https://github.com/[your-username]/[your-repo].git
cd [your-repo]
```

2. Place your Xenium data folders under `data/` and create a `results/` folder

3. Pre-download SingleR reference files and save to `data/`:

```r
saveRDS(celldex::HumanPrimaryCellAtlasData(), "data/hpca_ref.rds")
saveRDS(celldex::BlueprintEncodeData(),        "data/blueE_ref.rds")
```

4. In `load_libraries.R`, update the path lines:

```r
dataDir    <- "data/"
resultsDir <- "results/"
```

5. Set `use_toy_data <- FALSE` at the top of `analysis.R`

6. Run scripts in order:

```r
source("code/load_libraries.R")
source("code/functions.R")
source("code/analysis.R")
# Run figure_script_1.R twice — once per condition
Rscript code/figure_script_1.R Low
Rscript code/figure_script_1.R High
source("code/figure_script_2.R")
```

---

## Key design decisions

- **Batch integration:** Harmony is used to integrate 8 Xenium samples, run on PCA embeddings from SCTransform-normalised counts.
- **Cell type annotation:** Initial broad annotation uses SingleR with HumanPrimaryCellAtlas and BlueprintEncode references, followed by manual cluster renaming based on canonical marker genes. On toy data, UCell module scoring is used as an approximation.
- **CellChat:** Spatial CellChat is run separately per condition (VSTM5 Low/High, IFIT1B Low/High) using cell-type centroids and physical distance constraints (`interaction.range = 250`, `scale.distance = 7.4`).
- **Reproducibility:** All GitHub packages are pinned to specific commits. The Code Ocean capsule locks the full software environment via Docker.

---

## Data availability

Raw Xenium spatial transcriptomics data for ovarian cancer patient samples will be deposited in a public repository upon manuscript acceptance. Data access details will be updated here.

For enquiries regarding data access prior to publication, contact: [corresponding author email]

---

## Citation

If you use this code, please cite:

[Author list placeholder]. [Manuscript title placeholder]. *Nature Communications* (submitted).  
Code Ocean DOI: [PLACEHOLDER]  
GitHub: https://github.com/[your-username]/[your-repo]

---

## License

This code is released under the MIT License. See `LICENSE` for details.

---