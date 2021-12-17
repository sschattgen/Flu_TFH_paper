source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
library(parallel)
library(GSVA)
library(Matrix)
library(cerebroApp)
library(pheatmap)
library(ggplotify)
library(jcolors)
library(clusterProfiler)



# split by donor and recluster ====
Tfh_obj <- readRDS(Tfh_lineages_path)

donors <- unique(Tfh_obj@meta.data$donor)
years <- unique(Tfh_obj@meta.data$year)


donor_objs <- list()
for (i in donors){
  
  donor_obj_x <- subset(Tfh_obj, subset = donor == i)
  
  for (y in years){
    donor_obj <- subset(donor_obj_x, subset = year == y)
    donor_obj <- SeuratReprocess(donor_obj)
    donor_obj@meta.data$ident <- donor_obj@active.ident
    tlevels <- names(TimePal2)[which(names(TimePal2) %in% donor_obj@meta.data$time )]
    donor_obj@meta.data$time <- factor(donor_obj@meta.data$time, levels = tlevels)
    donor_obj@meta.data$year <- factor(donor_obj@meta.data$year, levels = y)
    donor_objs <- list.append(donor_objs, donor_obj)
  }
  
}

# Cerebro for gsva scanning ====

#gene sets
hallmark <- c("Z:/ResearchHome/Groups/thomagrp/home/sschattg/TOOLBOX/GSEA_gene_Sets/Hallmark_Human_Symbol_with_Remapping_MSi.gmt")
C2 <- c("Z:/ResearchHome/Groups/thomagrp/home/sschattg/TOOLBOX/GSEA_gene_Sets/c2.cp.kegg.v7.2.symbols.gmt")
C7 <- c("Z:/ResearchHome/Groups/thomagrp/home/sschattg/TOOLBOX/GSEA_gene_Sets/C7_Immune_Tcell_Human_Symbol_with_Remapping_MSi.gmt")
groups <- c('ident', 'time')
gslist <- list(hallmark, C2, C7)
names(gslist) <- c('hallmark', 'C2', 'C7')

#loop through by donor and 
for (z in seq_along(donor_objs)){
  
  for (ii in seq_along(gslist)){
    
    #gsva
    donor_objs[[z]] <- performGeneSetEnrichmentAnalysis(
      object = donor_objs[[z]],
      GMT_file = gslist[[ii]],
      groups = groups,
      thresh_p_val = 0.05,
      thresh_q_val = 0.1,
      name = names(gslist)[ii]
    )
    
  }
  
}


#significant hits 
donor_obj_hits <- list()
for (z in seq_along(donor_objs)){
  obj_hits <- list()
  for (x in names(gslist)){
    for (zz in groups){
      if( donor_objs[[z]]@misc$enriched_pathways[[x]][[zz]] == "no_gene_sets_enriched"){
        next
      } else{
        obj_hits <- list.append(obj_hits, donor_objs[[z]]@misc$enriched_pathways[[x]][[zz]])
      }
    }
  }
  donor_obj_hits <- list.append(donor_obj_hits, obj_hits)
}


