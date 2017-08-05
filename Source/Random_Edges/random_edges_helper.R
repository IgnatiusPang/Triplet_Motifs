#########################################################

randomly_add_or_remove_edges_from_network <- function ( network_edge_list, proportion, column_a, column_b,   
                                                        directed = FALSE, 
                                                        simple_network = TRUE,
                                                        maximum_num_of_tries = 10) {

  network_edge_list <- clean_up_network (  network_edge_list, column_a, column_b, directed=directed) 

  if ( proportion > 1) {
    
    network_edge_list <- randomly_add_edges_to_network(  network_edge_list,  proportion  - 1, column_a, column_b,   
                                                         directed = directed, 
                                                         simple_network = simple_network,
                                                         maximum_num_of_tries = maximum_num_of_tries  ) 
  } else if ( proportion < 1) {
    
    network_edge_list <- randomly_remove_edges_from_network   (  network_edge_list, 1- proportion  ) 
  } 

  return (network_edge_list ) 
}



######
# Floor of proportion_to_remove * number of edges in the network is used to determine the number of edges to remove 
randomly_remove_edges_from_network  <- function (  network_edge_list,  proportion_to_remove ) {
  
  
  if ( proportion_to_remove >= 1) {
    
    stop ( "randomly_remove_edges_from_undirected_network: proportion to remove needs to be smaller than 1")
    
  }
  
  sampled_network <- sample_frac(network_edge_list, size=1-proportion_to_remove, replace = FALSE)
  
  return( sampled_network)  
}


randomly_add_edges_to_network <- function (  network_edge_list,  proportion_to_add, column_a, column_b,   
                                                        directed = FALSE, 
											 			simple_network = TRUE,
                                                        maximum_num_of_tries = 10 ) {
  
  edge_list <- dplyr::select( network_edge_list, one_of ( c(column_a, column_b) ))
  
  multiple_edges_to_remove <- get_multiple_edges  (edge_list, column_a, column_b, directed=directed, duplicates=FALSE  ) 
  multiple_edges_to_remove_with_duplicates <- get_multiple_edges  (edge_list, column_a, column_b, directed=directed, duplicates=TRUE  ) 
  self_edges_to_remove     <- get_self_loops ( edge_list, column_a, column_b  ) 
  
  edge_list <- dplyr::setdiff( edge_list, multiple_edges_to_remove) %>% # Remove all multiple edges
    		   dplyr::setdiff( self_edges_to_remove) # Remove all self-loops

  edge_list <- dplyr::union ( edge_list, multiple_edges_to_remove) 
  
  num_edges <-  count ( edge_list )[[1]]
  
  num_edges_to_add <- ceiling ( proportion_to_add * num_edges )
  
  column_a_nodes <- unique( as.vector ( t(  as.data.frame(  edge_list[,column_a] ) ) ) ) 
  
  column_b_nodes <- unique( as.vector ( t(  as.data.frame(  edge_list[,column_b] ) ) ) ) 
  
  all_nodes <- unique( column_a_nodes, column_b_nodes)
  
  total_num_nodes <- length ( all_nodes)
  
  maximum_num_edges <- total_num_nodes * (total_num_nodes - 1) /2
  
  if (proportion_to_add > (maximum_num_edges - num_edges)/num_edges ) {
    stop ( "randomly_add_edges_to_undirected_network: Value for parameter 'proportion_to_add' is too large.")
  }

  num_of_tries <- 0 
  count_edges_left_to_add <- num_edges_to_add
  

  while ( count_edges_left_to_add > 0 & num_of_tries < maximum_num_of_tries  ) {
    
    source_nodes_to_add  <- c() 
    target_nodes_to_add  <- c()
    
    if ( directed == FALSE ) {
      source_nodes_to_add <- sample (all_nodes, count_edges_left_to_add, replace=TRUE )
      target_nodes_to_add <- sample (all_nodes, count_edges_left_to_add, replace=TRUE )
      
    } else {
      
      source_nodes_to_add <- sample (column_a_nodes, count_edges_left_to_add, replace=TRUE )
      target_nodes_to_add <- sample (column_b_nodes, count_edges_left_to_add, replace=TRUE )
    }
     
    edges_to_add <- cbind(source_nodes_to_add, target_nodes_to_add  )
    
    colnames(edges_to_add ) <- c( column_a, column_b)
    
    updated_edge_list <- rbind( edge_list, edges_to_add) 
    
    multiple_edges_updated <- get_multiple_edges  (updated_edge_list, column_a, column_b, directed=directed  ) 
    self_edges_updated     <- get_self_loops ( updated_edge_list, column_a, column_b  ) 
    
    updated_edge_list <- dplyr::setdiff( updated_edge_list, multiple_edges_updated) %>% # Remove all multiple edges
      					 dplyr::setdiff( self_edges_updated) # Remove all self-loops
    
    edge_list <- updated_edge_list
    
    
    current_num_edges <- count ( edge_list )[[1]]
    
    count_edges_left_to_add <- num_edges_to_add - ( current_num_edges -  num_edges)
    num_of_tries <- num_of_tries + 1
    
    # print ( count_edges_left_to_add)
  }
  
  ## Add original mulitple edges and self-loops back into the network 
  if ( simple_network == FALSE) {
	  edge_list <- dplyr::setdiff( edge_list, multiple_edges_to_remove) 
	  
	  edge_list <- dplyr::union ( edge_list, multiple_edges_to_remove_with_duplicates)   %>%
	  		    	dplyr::union( self_edges_to_remove) # Replace all self-loops
  }
  
  return ( edge_list)
}



clean_up_network <- function (  network_edge_list, column_a, column_b, directed=FALSE) {
	
	edge_list <- dplyr::select( network_edge_list, one_of ( c(column_a, column_b) ))
	
	multiple_edges_to_remove <- get_multiple_edges  (edge_list, column_a, column_b, directed=directed, duplicates=FALSE  ) 
	self_edges_to_remove     <- get_self_loops ( edge_list, column_a, column_b  ) 
	
	edge_list <- dplyr::setdiff( edge_list, multiple_edges_to_remove) %>% # Remove all multiple edges
		         dplyr::setdiff( self_edges_to_remove) # Remove all self-loops
	
	edge_list <- dplyr::union ( edge_list, multiple_edges_to_remove) 


	return(edge_list)	
}



#########################################################

