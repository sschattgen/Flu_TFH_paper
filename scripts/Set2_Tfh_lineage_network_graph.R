# Run Hyunjin's pipeline
source('Z:/ResearchHome/Groups/thomagrp/home/sschattg/bioinformatics_projects/Ali_Tfh/scripts/tfh_pkgs_paths_vars.R')
setwd(tfh_working_dir)
map <- purrr::map
# make clone lineage pivot table spread by time

make_wide_lineages <- function(df){
  df %>%
    select( clone_id, time_point, tissue, n  ) %>%
    mutate(time = gsub('_d','-d', time_point)) %>% 
    pivot_wider(values_from = 'n', names_from = 'time',  values_fill = 0) %>%
    mutate(across(where(is.integer), ~ nafill(., fill = 0))) %>%
    dplyr::rename(cell_type = tissue)
}

#lin network function
draw_lineage_network <- function(lineage_table,
                                 by.time=TRUE) {
  
  ### load library
  if(!require(igraph, quietly = TRUE)) {
    install.packages("igraph")
    require(igraph, quietly = TRUE)
  }
  if(!require(RedeR, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("RedeR")
    require(RedeR, quietly = TRUE)
  }
  
  ### get time points
  time_points <- colnames(lineage_table)[4:(ncol(lineage_table))]
  
  ### if by.time = TRUE
  if(by.time) {
    node_w_df <- lineage_table %>%
      group_by(time_point, cell_type) %>%
      tally() %>%
      mutate(name = paste(cell_type, gsub('_','-',time_point ), sep = '_'))
    node_weights <- node_w_df[['n']]
    names(node_weights) <- node_w_df[['name']]
    ### set node names
    node_names <- c(paste("PBMC", time_points, sep = "_"),
                    paste("FNA", time_points, sep = "_"))
    
    ### make adjacency matrix for network
    ### rows: outbound
    ### columns: inbound
    adj_mat <- matrix(0, length(node_names), length(node_names))
    rownames(adj_mat) <- node_names
    colnames(adj_mat) <- node_names
    
    for(r in rownames(adj_mat)) {
      for(c in colnames(adj_mat)) {
        ### outbound
        temp <- strsplit(r, split = "_", fixed = TRUE)[[1]]
        outbound_tissue <- temp[1]
        outbound_day <- temp[2]
        
        ### inbound
        temp <- strsplit(c, split = "_", fixed = TRUE)[[1]]
        inbound_tissue <- temp[1]
        inbound_day <- temp[2]
        
        ### tissue indicies
        ridx <- which(lineage_table$cell_type == outbound_tissue)
        cidx <- which(lineage_table$cell_type == inbound_tissue)
        
        ### fill out the table
        for(i in ridx) {
          for(j in cidx) {
            if(lineage_table$clone_id[i] == lineage_table$clone_id[j]) {
              if(lineage_table[i,outbound_day] > 0 && lineage_table[j,inbound_day] > 0) {
                adj_mat[r,c] <- adj_mat[r,c] + pull(lineage_table[j,inbound_day])
              }
              break;
            }
          }
        }
        
        ### print progress
        writeLines(paste(r, c))
      }
    }
    
    ### diagonal <- 0
    diag(adj_mat) <- 0
    
    ### remove nodes that do not have any connections (edges)
    remove_idx <- NULL
    for(i in 1:nrow(adj_mat)) {
      if((sum(adj_mat[i,]) == 0) && (sum(adj_mat[,i]) == 0)) {
        remove_idx <- c(remove_idx, i)
      }
    }
    adj_mat <- adj_mat[-remove_idx, -remove_idx]
    
    ### make an igraph
    g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
    coords <- layout_(g, as_tree())
    plot(g, layout = coords, vertex.label.size = 8)
    
    ### node and edge width + colors
    E(g)$width <- E(g)$weight
    E(g)$edgeColor <- "gray"
    V(g)$nodeSize <- sapply(V(g)$name, function(x) {
      return(max(adj_mat[,x]))
    })
    V(g)$nodeSize <- sapply(V(g)$name, function(x) {
      return(node_weights[[x]])
    })
    V(g)$nodeSize <- (V(g)$nodeSize / max(V(g)$nodeSize, na.rm = TRUE)) * 100
    V(g)$color <- sapply(V(g)$name, function(x) {
      if(grepl("PBMC", x, fixed = TRUE)) {
        return("#EA3426")
      } else if(grepl("FNA", x, fixed = TRUE)) {
        return("#0C1BFF")
      } else {
        return("white")
      }
    })
    
    ### edge sizes should be smaller than node sizes
    if(max(E(g)$width) > (0.3 * max(V(g)$nodeSize))) {
      E(g)$width <- (E(g)$width / max(E(g)$width, na.rm = TRUE)) * (0.3 * max(V(g)$nodeSize)) 
    }
    
    ### increase the vertex label size
    V(g)$label.cex <- 2
    
    ### load RedeR screen and plot the graph
    rdp<-RedPort()
    calld(rdp)
    addGraph(rdp,g, layout.kamada.kawai(g))
    
    ### add legends
    # color
    addLegend.color(rdp, colvec=c("#EA3426", "#0C1BFF"), labvec=c("PBMC", "FNA"), title="",
                    vertical=TRUE, position="bottomleft", dxborder=50, dyborder=400, size=40, ftsize=30)
    
    # size
    circleLabel <- floor(seq(min(V(g)$nodeSize),max(V(g)$nodeSize),(max(V(g)$nodeSize) - min(V(g)$nodeSize))/4))
    circleSize <- (circleLabel / max(circleLabel)) * 100
    diag(adj_mat) <- NA
    circleLabel <- floor(seq(min(adj_mat, na.rm = TRUE), max(adj_mat, na.rm = TRUE),
                             ((max(adj_mat, na.rm = TRUE) - min(adj_mat, na.rm = TRUE))/4)))
    ### circle size in the legend should be at least 1
    if(circleSize[1] == 0) {
      circleLabel[1] <- 1
      circleSize[1] <- (1 / max(circleLabel)) * 100
    }
    addLegend.size(rdp,sizevec=circleSize,labvec=circleLabel,title="Clone Size",
                   position="bottomleft", dxborder=50, dyborder=100, ftsize=24)
    
  } 
  
}
draw_lineage_network2 <- function(lineage_table,
                                 by.time=TRUE) {
  
  ### load library
  if(!require(igraph, quietly = TRUE)) {
    install.packages("igraph")
    require(igraph, quietly = TRUE)
  }
  if(!require(RedeR, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("RedeR")
    require(RedeR, quietly = TRUE)
  }
  
  ### get time points
  time_points <- colnames(lineage_table)[4:(ncol(lineage_table))]
  
  ### if by.time = TRUE
  if(by.time) {
    
    ### set node names
    node_names <- c(paste("PBMC", time_points, sep = "_"),
                    paste("FNA", time_points, sep = "_"))
    
    ### make adjacency matrix for network
    ### rows: outbound
    ### columns: inbound
    adj_mat <- matrix(0, length(node_names), length(node_names))
    rownames(adj_mat) <- node_names
    colnames(adj_mat) <- node_names
    
    for(r in rownames(adj_mat)) {
      for(c in colnames(adj_mat)) {
        ### outbound
        temp <- strsplit(r, split = "_", fixed = TRUE)[[1]]
        outbound_tissue <- temp[1]
        outbound_day <- temp[2]
        
        ### inbound
        temp <- strsplit(c, split = "_", fixed = TRUE)[[1]]
        inbound_tissue <- temp[1]
        inbound_day <- temp[2]
        
        ### tissue indicies
        ridx <- which(lineage_table$cell_type == outbound_tissue)
        cidx <- which(lineage_table$cell_type == inbound_tissue)
        
        ### fill out the table
        for(i in ridx) {
          for(j in cidx) {
            if(lineage_table$clone_id[i] == lineage_table$clone_id[j]) {
              if(lineage_table[i,outbound_day] > 0 && lineage_table[j,inbound_day] > 0) {
                adj_mat[r,c] <- adj_mat[r,c] + pull(lineage_table[j,inbound_day])
              }
              break;
            }
          }
        }
        
        ### print progress
        writeLines(paste(r, c))
      }
    }
    
    ### diagonal <- 0
    diag(adj_mat) <- 0
    
    ### remove nodes that do not have any connections (edges)
    remove_idx <- NULL
    for(i in 1:nrow(adj_mat)) {
      if((sum(adj_mat[i,]) == 0) && (sum(adj_mat[,i]) == 0)) {
        remove_idx <- c(remove_idx, i)
      }
    }
    adj_mat <- adj_mat[-remove_idx, -remove_idx]
  }
  return(adj_mat)
}
#make lineage input ====
Clonedf <- read.delim(clone_df_path, stringsAsFactors = F) %>%
  mutate(Tfh_type2 = ifelse(Tfh_type %notin% c("GC","pre/memory","IL10 TFH", "high ISG","cycling") ,NA, 'Tfh' )) %>%
  mutate(time_point = factor(time_point, level=names(TimePal2))) %>%
  filter(!is.na(day))

tfh_freq_allivium_df <- Clonedf %>%
  select(donor, time_point, tissue, clone_id, Tfh_type2) %>%
  group_by_all() %>%
  tally() %>%
  ungroup(clone_id, Tfh_type2) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  filter(!is.na(Tfh_type2)) %>%
  mutate(time_point = factor( time_point, levels = names(TimePal2))) 

split_clone_freq_list <- tfh_freq_allivium_df %>%
  group_by(donor) %>%
  group_split()

lineage_dfs <- map( split_clone_freq_list, ~make_wide_lineages(.))

# draw lineages====
draw_lineage_network(lineage_dfs[[1]])
draw_lineage_network(lineage_dfs[[2]])
draw_lineage_network(lineage_dfs[[3]])
draw_lineage_network(lineage_dfs[[4]])




#mat <- draw_lineage_network2(lineage_dfs[[1]])
#write_tsv(as.data.frame(mat), '~/d4_adj_matrix.adj')