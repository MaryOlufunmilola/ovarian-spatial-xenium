#!/bin/bash
set -e  # stop immediately if any script throws an error
 
# ----------------------------------------------------------
# run.sh — Code Ocean Capsule Entry Point
# Scripts execute in this order every Reproducible Run.
# ----------------------------------------------------------
 
# Step 1: Load libraries and set paths
Rscript /code/load_libraries.R
 
# Step 2: Load shared functions
Rscript /code/functions.R
 
# Step 3: Main analysis — processing, integration, annotation
#         Saves /results/Seurat_object.rds
Rscript /code/analysis.R
 
# Step 4: CellChat pipeline — one condition at a time to manage RAM (>64 GB each)
Rscript /code/figure_script_1.R Low
Rscript /code/figure_script_1.R High
 
# Step 5: Manuscript figures (requires both CellChat RDS files from Step 4)
Rscript /code/figure_script_2.R