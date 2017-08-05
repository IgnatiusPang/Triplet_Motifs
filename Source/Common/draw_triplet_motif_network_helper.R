
# Function: draw_triplet_motif_network

## The input has to following the following format: 

#   oln_id_a oln_id_b oln_id_c type_ac type_bc genetic_interaction_score
# 1  YBL007C  YAL023C  YKL112W      tu      tu                   -0.1390
# 2  YBL023C  YAR007C  YDL056W      tu      tu                   -0.1556

draw_triplet_motif_network <- function ( triplet_motifs_costanzo, network_name ) { 
	
	#######################################################################################################################
	
	## Add the list of nodes 
	list_of_nodes <- unique( c(triplet_motifs_costanzo[,1], triplet_motifs_costanzo[,2], triplet_motifs_costanzo[,3]))
	
	graph_nel_object <- new ('graphNEL', edgemode='undirected')
	
	graph_nel_object <- graph::addNode( as.character( list_of_nodes) , graph_nel_object)
	
	#######################################################################################################################
	
	# Create Cytoscape window
	cw <- CytoscapeWindow ( network_name, graph=graph_nel_object, overwrite=TRUE)
	graph_nel_object <- cw@graph
	
	#######################################################################################################################
	
	## Add edges
	
	graph_nel_object <- initEdgeAttribute (graph_nel_object, "edgeType", "char", "gi")
	
	gene_a <- as.character( as.vector(c( t(triplet_motifs_costanzo[,1]) ,
										 t(triplet_motifs_costanzo[,1])    , 
										 t(triplet_motifs_costanzo[,2]) ) ))
	
	gene_b <- as.character( as.vector(c( t(triplet_motifs_costanzo[,2]) ,
										 t(triplet_motifs_costanzo[,3])    , 
										 t(triplet_motifs_costanzo[,3]) ) ))
	
	graph_nel_object <- graph::addEdge(gene_a, gene_b, graph_nel_object) 
	
	edge_type_attributes <- as.character( as.vector( c( 
		rep( "gi", length(triplet_motifs_costanzo[,1] )),
		t(triplet_motifs_costanzo[, c("type_ac")]),
		t(triplet_motifs_costanzo[, c("type_bc")])
	)  ) ) 
	
	edgeData( graph_nel_object, gene_a, gene_b, "edgeType") <- edge_type_attributes
	
	#######################################################################################################################
	
	### Create the network visualization
	
	cw <- setGraph(cw, graph_nel_object)
	
	displayGraph(cw)
	
	layoutNetwork(cw, layout.name = "force-directed")
	
	return (cw )
}

#######################################################################################################################
#######################################################################################################################

# > getArrowShapes(cw)
# [1] "DIAMOND_SHORT_2" "DELTA_SHORT_2"   "HALF_BOTTOM"     "ARROW_SHORT"     "T"               "DELTA_SHORT_1"   "HALF_TOP"        "DELTA"           "DIAMOND"        
# [10] "CIRCLE"          "DIAMOND_SHORT_1" "ARROW"           "NONE" 
blue <- '#0000FF'
red <- '#FF0000'
black <- '#000000'
green <- '#00FF00'
grey <- '#808080'
	
triplet_motif_network_display_style <- function(cw) {

	setEdgeColorRule(cw, 'edgeType',
					 control.points=c( 'ku', 'kd', 'gi', 'tu', 'td', 'p'), 
					 colors=c(blue, blue, green, red, red, black), 
					 mode='lookup', default.color=grey)
	
	setEdgeSourceArrowRule(cw, 'edgeType',
						   attribute.values=c( 'ku', 'kd', 'gi', 'tu', 'td', 'p'), 
						   arrows=c('ARROW', 'NONE', 'NONE', 'ARROW', 'NONE', 'NONE'), 
						   default='NONE')
	
	setEdgeTargetArrowRule(cw, 'edgeType',
						   attribute.values=c( 'ku', 'kd', 'gi', 'tu', 'td', 'p'), 
						   arrows=c('NONE', 'ARROW', 'NONE', 'NONE', 'ARROW', 'NONE'),  
						   default='NONE')
	
	setEdgeSourceArrowColorRule(cw, 'edgeType',
								control.points=c( 'ku', 'kd', 'gi', 'tu', 'td', 'p'), 
								colors=c('#0000FF', blue, green, red, red, black), 
								mode='lookup', default.color=grey)
	
	setEdgeTargetArrowColorRule(cw, 'edgeType',
								control.points=c( 'ku', 'kd', 'gi', 'tu', 'td', 'p'), 
								colors=c(blue, blue, green, red, red, black),  
								mode='lookup', default.color=grey)
	
	if ( 'gene_name' %in% noa.names(getGraph(cw))) {
		
		setNodeLabelRule(cw, 'gene_name')
		
	}
	

	if ( 'Domain..predominant.' %in% noa.names(getGraph(cw))) {
		
		setNodeColorRule(cw, "Domain..predominant.", control.points=1:17, colors= c( "#CCCCCC"
            , "#000000"																					 
			,"#00FFCC"	
			,"#000099"	
			,"#CC9900"	
			,"#990099"	
			,"#999900"	
			,"#00FFCC"	
			,"#CC33FF"	
			,"#00CC99"	
			,"#33CC00"	
			,"#996600"	
			,"#66FF99"	
			,"#CC0033"	
			,"#6666FF"	
			,"#33FF99"	
			,"#CCCCFF"), mode='lookup')
	}
	
}


#######################################################################################################################

draw_pairwise_interactions <- function ( edge_list, network_name ) { 
	
	#######################################################################################################################
	## Add the list of nodes 
	
	list_of_nodes <- unique( c(edge_list[,1], edge_list[,2]))
	
	graph_nel_object <- new ('graphNEL', edgemode='undirected')
	
	graph_nel_object <- graph::addNode( as.character( list_of_nodes) , graph_nel_object)
	
	#######################################################################################################################
	# Create Cytoscape window
	
	cw <- CytoscapeWindow ( network_name, graph=graph_nel_object, overwrite=TRUE)
	graph_nel_object <- cw@graph
	
	#######################################################################################################################
	## Add edges
	
	gene_a <- as.character( as.vector(c( t(edge_list[,1]) ) ))
	
	gene_b <- as.character( as.vector(c( t(edge_list[,2]) ) ))
	
	if ( length ( edge_list[1,]) > 2 ) {
		
		
		for ( i in 3:length(edge_list[1,] ) ) { 
			
			edge_attribute_name <- colnames(edge_list)[i]
			
			if( edge_attribute_name == 'edgetype') {
				edge_attribute_name <- 'edgeType'
				
			}
			
			graph_nel_object <- initEdgeAttribute (graph_nel_object, edge_attribute_name, "char", edge_list[1,i])
		
			graph_nel_object <- graph::addEdge(gene_a, gene_b, graph_nel_object) 
			
			edge_type_attributes <- as.character( as.vector( t( edge_list[,i] )  ) ) 
			
			edgeData( graph_nel_object, gene_a, gene_b, edge_attribute_name) <- edge_type_attributes
		
		}
		
	}

	#######################################################################################################################
	
	### Create the network visualization
	
	cw <- setGraph(cw, graph_nel_object)

	
	return (cw )
}


#######################################################################################################################

get_attribute_type_and_default_value <- function( r_type ) {
	
	attribute.type <- ""
	default.value <- ""
	
	if ( r_type %in% c('character')) {
		attribute.type <- "char"
		default.value <- "undefined"
	} else if ( r_type %in% c( 'integer', 'int')) {
		attribute.type <- "integer"
		default.value <- NA
	} else if ( r_type %in% c('double', 'numeric')) {
		attribute.type <- "numeric"
		default.value <- NA
	}
	
	return ( list(attribute.type= attribute.type, default.value=default.value ))
}

## Add node attributes for existing nodes
## cw - cytoscape workspace
## node_attributes - table containing the node attributes
## node_key_column - the primary key for the nodes

add_node_attributes <- function ( cw, node_attributes, node_key_column ) {
	
	graph_nel_object <- cw@graph
	
	existing_nodes <- nodes(graph_nel_object)
	
	## Add attributes for nodes that are already present in the network 
	
	nodes_to_add_attributes <- intersect (  existing_nodes, t( node_attributes[, node_key_column] )  ) %>% as.data.frame()
	
	colnames( nodes_to_add_attributes) <- node_key_column
	
	# print ( existing_nodes)
	
	# Subset of the table 
	subset_attribute_table <- inner_join( node_attributes, nodes_to_add_attributes, by=node_key_column)
	
	# Update the order of nodes_to_add_attributes	
	nodes_to_add_attributes <- as.vector ( t(subset_attribute_table[,node_key_column]) )
	
	
	# print ( head( subset_attribute_table) ) 
	
	# For each column 
	for(i in seq_along(subset_attribute_table[1,])) {
		
		current_column_name   <- colnames(subset_attribute_table)[i]
		typeof_current_column <- typeof( subset_attribute_table[1,i])
		type_and_default_value <-  get_attribute_type_and_default_value ( typeof_current_column )
		
		attribute.type 	      <- type_and_default_value$attribute.type
		default.value         <- type_and_default_value$default.value

		if ( current_column_name != node_key_column ) {
			
			graph_nel_object <- initNodeAttribute( graph=graph_nel_object, 
												   attribute.name=current_column_name, 
												   attribute.type=attribute.type, 
												   default.value=default.value )
			
			the_attributes_to_add <- as.vector ( t(subset_attribute_table[,i]) ) 
			
			nodeData (graph_nel_object, nodes_to_add_attributes, current_column_name) <- the_attributes_to_add
			
			# print ( nodes_to_add_attributes)
			# print( the_attributes_to_add )
		}
		
	}
	
	cw <- setGraph(cw, graph_nel_object)

	
	return(cw)
	
}

#######################################################################################################################

edge_list_network_display_style <- function(cw) {
	
	setNodeLabelRule(cw, 'gene_name')
	
	setDefaultNodeColor(cw, '#66CCFF')
	
	setNodeColorRule(cw, 'is_essential', control.points=as.integer(c(0,1)), colors=c('#66CCFF', '#FF6666'), mode='lookup', default.color=grey)
	
	setEdgeColorRule(cw, 'edgeType',
					 control.points=c( 'ppi_n_ngi'), 
					 colors=c(black), 
					 mode='lookup', default.color=grey)
}

#######################################################################################################################