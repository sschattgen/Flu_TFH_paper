
if(!require(data.table, quietly = TRUE)) {
  install.packages(data.table, repos = 'https://cloud.r-project.org/')
} 

if(!require(tidyverse, quietly = TRUE)) {
  install.packages(tidyverse, repos = 'https://cloud.r-project.org/')
} 

if(!require(Seurat, quietly = TRUE)) {
  install.packages(Seurat, repos = 'https://cloud.r-project.org/')
} 

if(!require(patchwork, quietly = TRUE)) {
  install.packages(patchwork, repos = 'https://cloud.r-project.org/')
} 

if(!require(jcolors, quietly = TRUE)) {
  install.packages(jcolors, repos = 'https://cloud.r-project.org/')
} 

if(!require(viridis, quietly = TRUE)) {
  install.packages(viridis, repos = 'https://cloud.r-project.org/')
} 

if(!require(ggthemes, quietly = TRUE)) {
  install.packages(ggthemes, repos = 'https://cloud.r-project.org/')
} 

if(!require(ggpubr, quietly = TRUE)) {
  install.packages(ggpubr, repos = 'https://cloud.r-project.org/')
} 
if(!require(rlist, quietly = TRUE)) {
  install.packages(rlist, repos = 'https://cloud.r-project.org/')
} 
if(!require(ggrepel, quietly = TRUE)) {
  install.packages(ggrepel, repos = 'https://cloud.r-project.org/')
}
if(!require(ggplotify, quietly = TRUE)) {
  install.packages(ggplotify, repos = 'https://cloud.r-project.org/')
} 
if(!require(wesanderson, quietly = TRUE)) {
  install.packages(wesanderson, repos = 'https://cloud.r-project.org/')
} 
if(!require(ggalluvial, quietly = TRUE)) {
  install.packages(ggalluvial, repos = 'https://cloud.r-project.org/')
} 
if(!require(pals, quietly = TRUE)) {
  install.packages(pals, repos = 'https://cloud.r-project.org/')
} 
if(!require(igraph, quietly = TRUE)) {
  install.packages(igraph, repos = 'https://cloud.r-project.org/')
} 
if(!require(ggplotify, quietly = TRUE)) {
  install.packages(ggplotify, repos = 'https://cloud.r-project.org/')
} if(!require(pheatmap, quietly = TRUE)) {
  install.packages(pheatmap, repos = 'https://cloud.r-project.org/')
} 
if(!require(Peptides, quietly = TRUE)) {
  install.packages(Peptides, repos = 'https://cloud.r-project.org/')
} 
if(!require(RedeR, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("RedeR")
}

if(!require(cerebroApp, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install('romanhaa/cerebroApp')
}
if(!require(org.Hs.eg.db, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}
if(!require(GSVA, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install('GSVA')
}