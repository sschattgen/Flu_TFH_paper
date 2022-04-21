#!/bin/bash
cd path/to/repository
mkdir outs
Rscript scripts/setup.R

# Preprocessing
#python scripts/TwoYear_VDJ_parsing.py
#Rscript scripts/TwoYear_gen_seurat.R

# T cell analysis
Rscript scripts/TwoYear_Tcell_annotation.R
Rscript scripts/TwoYear_10x_clones_annotate.R
Rscript scripts/TwoYear_Tcell_analysis.R
Rscript scripts/TwoYear_Tfh_subset.R
Rscript scripts/TwoYear_Tfh_clonotypes.R
Rscript scripts/TwoYear_Tfh_phenotypes.R 
Rscript scripts/TwoYear_Tcell_metabolism.R
Rscript scripts/TwoYear_Tfh_time.R
Rscript scripts/TwoYear_Donor_321-05_pickedClones.R 
Rscript scripts/TwoYear_Donor_321-05_pickedClones_gex_variation.R
Rscript scripts/TwoYear_Donor_321-05_pickedClones_modules.R
Rscript scripts/TwoYear_B_vs_T_corr.R
Rscript scripts/TwoYear_run_HK_pipeline.R
