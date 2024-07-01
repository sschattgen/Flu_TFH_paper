# TwoYear Tfh phenotypes over time

source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)

rename <- dplyr::rename
# import T cells Seurat object and df ----
Tfh_obj <- readRDS(Set2_integrated_Tfh_path)
gene_set <- read.delim( './10x/outs/mart_export.txt') %>%
  pull(Human.gene.name) %>%
  unique()

TFH2 <- subset(Tfh_obj, Tfh_type != 'Treg') 




TFH2 <- AddModuleScore(TFH2, list(full_mod = gene_set), nbin = 15)

vln2 <- VlnPlot(TFH2,'Cluster1', group.by = 'Tfh_type',pt.size =0) +
  scale_fill_manual(values = TfhPal) +
  xlab(NULL) +
  labs(subtitle = 'Full vs. Prog TFH Module Score\nPodestà et al., 2023', title = NULL) +
  NoLegend() 
ggsave('./10x/outs/Set2_podesta_mod_score_vln.png', 
       plot = vln2, width = 4, height = 3)
ggsave('./10x/outs/Set2_podesta_mod_score_vln.pdf', 
       plot = vln2, width = 4, height = 3)