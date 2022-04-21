# Run Hyunjin's pipeline
source('scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
source('scripts/TFH_Analyses_TwoYear.R')
tfh_analyses_both_donors(Seurat_RObj_path = Tfh_lineages_path, 
                         outputDir = 'outs/')

