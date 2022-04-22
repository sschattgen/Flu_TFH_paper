#!/bin/bash
mkdir outs
Rscript scripts/setup.R

# Preprocessing
#python scripts/TwoYear_VDJ_parsing.py
#Rscript scripts/TwoYear_gen_seurat.R

## T cell analysis ##

#Rebuild the T cell Seurat object and annotate the clusters
Rscript scripts/TwoYear_Tcell_annotation.R
# Match 10X clones to the bulk and single-cell TCR data
Rscript scripts/TwoYear_10x_clones_annotate.R
# Basic stats and info shown in figure 2.
Rscript scripts/TwoYear_Tcell_analysis.R
# Subsets on the Tfh cells and makes some plots for figure 3.
Rscript scripts/TwoYear_Tfh_subset.R
# clonality and diversity stats from figure 6.
Rscript scripts/TwoYear_Tfh_clonotypes.R
# DEG id for TFH subsets for figure 3.
Rscript scripts/TwoYear_Tfh_phenotypes.R
# Gene set enrichment analysis for figure 4. 
Rscript scripts/TwoYear_Tcell_metabolism.R
# time plot in figure 5 and other unused things. 
Rscript scripts/TwoYear_Tfh_time.R
# identification of the picked clones for figure 7. 
Rscript scripts/TwoYear_Donor_321-05_pickedClones.R 
# picked clones gex analysis for figure 7. 
Rscript scripts/TwoYear_Donor_321-05_pickedClones_gex_variation.R
# picked clones gex analysis for figure 7.
Rscript scripts/TwoYear_Donor_321-05_pickedClones_modules.R
# figure 3E.
Rscript scripts/TwoYear_B_vs_T_corr.R
# psudotime and network graph of clones.
Rscript scripts/TwoYear_run_HK_pipeline.R
