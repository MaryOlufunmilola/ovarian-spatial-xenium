# Xenium Spatial Transcriptomics — Ovarian Cancer

R code for the Xenium spatial transcriptomics component of a multi-modal ovarian cancer study.  
**Authors:** Funmi Oyebamiji (analysis), Eleanor Paskus (figures)  
**Journal:** Nature Communications (submitted)

---

## Repository contents

| File | Description |
|---|---|
| `code/load_libraries.R` | Package installation and path setup |
| `code/functions.R` | Shared utility functions |
| `code/analysis.R` | Main analysis — processing, integration, annotation |
| `environment/postInstall` | GitHub package installer |
| `run.sh` | Execution entry point |

> Manuscript figure scripts (`figure_script_1.R`, `figure_script_2.R`) by Eleanor Paskus are available upon request.

---

## Requirements

R ≥ 4.4.3. All packages installed via `load_libraries.R`.  
GitHub-only packages (`ComplexHeatmap`, `circlize`, `CellChat`, `presto`) installed via `environment/postInstall`.

> ⚠️ **Computational requirements:** The CellChat spatial pipeline requires >64 GB RAM and several hours of runtime per condition. The full analysis is not suitable for standard cloud compute environments.

---

## Running the analysis

Update `dataDir` and `resultsDir` in `load_libraries.R` to point to your local data and results folders, then run scripts in order:

```r
source("code/load_libraries.R")
source("code/functions.R")
source("code/analysis.R")
```

Figure scripts require CellChat objects saved from `analysis.R` and are available upon request.

---

## Data availability

Raw Xenium data are not included due to patient privacy. Data available upon reasonable request to the corresponding author subject to IRB approval.

---

## License

[MIT License](https://opensource.org/licenses/MIT)
