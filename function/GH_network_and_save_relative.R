##### FONCTION TO PRODUCE A CO-OCCURRENCE NETWORK AND SAVE A ASV FILE #####

# seq = dataframe of the relative abundances (rows = ASV, columns = samplings)
# rho = Pearson's rho correlation coefficient threshold
# pvalue = p-value threshold for the significativity of the correlation

network_and_save_relative <- function(seq = seq,
                             rho = rho,
                             pvalue = p) {
  ## Relative abundances 
  seq_names <- seq$Representative_Sequence
  seq <- as.data.frame(lapply(seq[,-c(1,length(seq)-1,length(seq))], function(x) if(is.character(x)) as.numeric(x) else x))
  seq_r <- make_relative(as.matrix(t(seq)))
  seq_rh <- decostand(seq_r, method = "hellinger") 
  seq_rh <- as.data.frame(seq_rh)
  base::colnames(seq_rh) <- seq_names
  
  binary <- seq_rh[, colSums(seq_rh) != 0] # remove ASV in 0 sampling (not case here)
  
  ##  Correlation matrix
  corr_binary <- rcorr(as.matrix(binary), type = 'spearman')$r # Spearman's correlation coefficient
  corr_binary[corr_binary <= rho] <- 0
  
  pval_binary <- rcorr(as.matrix(binary), type = 'spearman')$P # Associated p-value
  
  for (i in 1:dim(corr_binary)[1]) {
    for (j in 1:dim(corr_binary)[1]) {
      corr_binary[i,j] <- ifelse(pval_binary[i,j] < pvalue, corr_binary[i,j], 0)
    }
  }
  remove(i, j)
  
  ## Create the network
  g1 <- graph_from_adjacency_matrix(corr_binary, mode = 'undirected', weighted = T, diag = FALSE)

  ## Network parameters
  # Degree 
  deg_node <- as.data.frame(igraph::degree(g1)) # All nodes degrees
  deg_node <- rownames_to_column(deg_node)
  colnames(deg_node) <- c('Representative_Sequence', 'Degree')
  
  deg_network <- mean(deg_node$Degree) # Mean of all the nodes degrees
  
  # Betweenness
  btw_node <- as.data.frame(betweenness(g1, directed = FALSE)) # All the nodes betweenness
  btw_node <- rownames_to_column(btw_node)
  colnames(btw_node) <- c('Representative_Sequence', 'Betweenness')

  btw_network <- mean(btw_node$Betweenness) # Mean of all the nodes betweenness

  # Centrality = Closeness
  cln_node <- as.data.frame(closeness(g1)) # All the nodes centralities
  cln_node <- rownames_to_column(cln_node)
  colnames(cln_node) <- c('Representative_Sequence', 'Centrality')

  cln_network <- mean(cln_node$Centrality) # Mean of all the nodes centralities

  # Transitivity = Clustering coefficient
  trs_node <- as.data.frame(transitivity(g1, type = 'local')) # All nodes transitivities
  
  trs_network <- transitivity(g1, type = 'global') # Mean of all the nodes transitivities
  
  # Other parameters
  dst_network <- edge_density(g1) # ratio # edges / # possible edges
  dia_network <- diameter(g1) # Diameter = length of the longest path
  ne_netork <- ecount(g1) # # of edges
  nv_network <- vcount(g1) # # of nodes
  spl_network <- mean_distance(g1) # Mean of the average path lengths between 2 nodes
  con_network <- vertex.connectivity(g1) # Connectivity = number of nodes to remove to isolate the node from the rest of the graph, here = 0 so graph not strongly connected
  
  ## Merging
  topo_network <- data.frame(ne_netork, nv_network, dst_network,
                             deg_network, dia_network, spl_network,
                             btw_network, cln_network, trs_network,
                             con_network)
  
  topo_node <- data.frame(deg_node, btw_node$Betweenness, cln_node$Centrality, trs_node)
  colnames(topo_node)[3:5] <- c('Betweenness', 'Centrality', 'Transitivity')
  
  results <- list(graph = g1, 
                  network = topo_network, 
                  node = topo_node)
  
  return(results)
}