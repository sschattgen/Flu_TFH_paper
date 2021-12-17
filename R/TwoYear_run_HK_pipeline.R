# Run Hyunjin's pipeline
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

source('./scripts/final_scripts/TFH_Analyses_TwoYear.R')

tfh_analyses_both_donors(Seurat_RObj_path = Tfh_lineages_path, 
                         outputDir = './10x/outs/HK/')

