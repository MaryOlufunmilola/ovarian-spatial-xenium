# Spatial Transcriptomics Analysis of Ovarian Cancer Patient Samples

> **Manuscript title:** [PLACEHOLDER — update when accepted]  
> **Journal:** Nature Communications (submitted)  
> **Code authors:** Funmi Oyebamiji, Eleanor Paskus

---

## Overview

This repository contains the R code by Funmi Oyebamiji for the **Xenium spatial transcriptomics** component of a multi-modal study of ovarian cancer patient samples.
This codebase covers spatial transcriptomics processing, clustering, cell-type annotation, and figure generation.

A fully reproducible computational capsule is available on Code Ocean:  
**Code Ocean DOI:** `[PLACEHOLDER — add after capsule is published]`

---

## Repository structure

```
.
├── code/
│   ├── load_libraries.R    # Package installation and path setup 
│   ├── functions.R         # Shared utility functions 
│   ├── analysis.R          # Main analysis — processing, integration, annotation 
├── environment/
│   └── postInstall         # GitHub-sourced R package installer (for Code Ocean)
├── .gitignore
├── run.sh                  # Code Ocean entry point — executes scripts in order
└── README.md
```

> **Note on figure scripts:** The manuscript figure scripts (`figure_script_1.R`, `figure_script_2.R`) were written by Eleanor Paskus, modified by Funmi Oyebamiji and are available in the Code Ocean capsule but are not hosted in this GitHub repository.

> **Data not included.** Raw and processed Xenium files are not hosted here due to size and patient privacy. See the Data Availability section below.

---

## Requirements

### R version
R ≥ 4.4.3 recommended

### CRAN / Bioconductor packages
All packages are loaded via `load_libraries.R`. Key dependencies include:

| Package | Source | Purpose |
|---|---|---|
| Seurat (v5) | CRAN | Spatial object handling, clustering |
| scCustomize | CRAN | Seurat visualisation helpers |
| SingleR | CRAN | Automated cell-type annotation |
| UCell | CRAN | Gene signature scoring |
| harmony | CRAN | Batch integration |
| CellChat | GitHub (jinworks) | Cell-cell communication |
| ComplexHeatmap | GitHub (jokergoo) | Heatmap figures |
| BiocParallel | Bioconductor | Parallel processing |
| org.Hs.eg.db | Bioconductor | Human gene annotation |

### GitHub-only packages
The following are installed via `environment/postInstall` (or manually with `devtools::install_github`):

```r
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jokergoo/circlize")
devtools::install_github("jinworks/CellChat")
devtools::install_github("immunogenomics/presto")
```

---

## Reproducing the analysis

### On Code Ocean (recommended)
1. Open the Code Ocean capsule at the DOI above
2. Xenium data files are attached to the capsule under `/data`
3. Click **Reproducible Run** — scripts execute automatically via `run.sh`

> **Important — two-run requirement for `figure_script_1.R`**
>
> The CellChat spatial pipeline (`figure_script_1.R`) must be run **twice**, once per condition (VSTM5 Low, VSTM5 High) due to RAM constraints. Each run requires >64 GB of RAM; running both conditions simultaneously would exceed available memory.
>
> The automated **Reproducible Run** (via `run.sh`) executes the Low condition by default. To reproduce the High condition:
>
> 1. In the Code Ocean capsule, open a **Terminal**
> 2. In the terminal, run with the condition passed as an argument — no file editing needed:
>    ```bash
>    Rscript /code/figure_script_1.R High
>    ```
> 4. This saves `/results/vstm5high.rds`, which `figure_script_2.R` depends on
>
> Both `/results/vstm5low.rds` and `/results/vstm5high.rds` must exist before running `figure_script_2.R`.

### Locally
1. Clone this repository:
   ```bash
   git clone https://github.com/[your-username]/[your-repo].git
   cd [your-repo]
   ```
2. Place your Xenium data files in a local `data/` folder and create a `results/` folder
3. In `load_libraries.R` change the two path lines to point to your local folders:
   ```r
   dataDir    <- "data/"
   resultsDir <- "results/"
   ```
4. Run scripts in order:
   ```r
   source("code/load_libraries.R")
   source("code/functions.R")
   source("code/analysis.R")
   ```
5. `figure_script_1.R` and `figure_script_2.R` (Ela Pask) are available in the Code Ocean capsule

---

## Data availability

Raw Xenium spatial transcriptomics data for ovarian cancer patient samples will be deposited in a public repository upon manuscript acceptance. Data access details will be updated here.

For enquiries regarding data access prior to publication, contact: [PLACEHOLDER — corresponding author email]

---

## Citation

If you use this code, please cite:

> [Author list placeholder]. *[Manuscript title placeholder]*. Nature Communications (submitted). Code Ocean DOI: `[PLACEHOLDER]`. GitHub: `https://github.com/[your-username]/[your-repo]`.

---

## License

This code is released under the [MIT License](https://opensource.org/licenses/MIT). See `LICENSE` for details.

---