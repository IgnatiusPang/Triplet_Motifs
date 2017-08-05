### Script: count_triplet_motifs_helper.R
### Author: Ignatius Pang 
### Date: 20-5-2016
### Description: Count the observed number of triplet motifs. Randomize the respective networks and count the number triplet motifs, repeat randomization 1000 times.
### Obtain the bootstrap p-value. This script contain the helper functions for the main script 'count_triplet_motifs_negative_interactions.R'.

## Project space:
# cd '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/'
# psql sbi_triplet_motifs 

#########################################################
## Insert directed network into both directions into the edge list table, so the direction that tables join together won't affect the analyses

## Inputs:
## input_table - must have 4 columns only, which are: gene_a, gene_b, interaction type, interaction type abbreviation
## gene_a - column name for input table that cotains gene A of the interaction
## gene_b - column name for input table that cotains gene B of the interaction
## interaction_type_forward - information on the type of interaction for A  ---> B 
## interaction_type_reverse - information on the type of interaction for B <---> A
## interaction_abbrev_forward - abbereviated information on the type of interaction for A ---> B 
## interaction_abberv_reverse - abbereviated information on the type of interaction for B ---> A 
## col_a - column name for gene A in the output table
## col_b - column name for gene B in the output table
## col_type - column name for the interaction type in the output table
## col_abbrev - column name for the interaction type abbreviation in the output table

## Output table:
## Output table has 4 columns, denoted by col_a, col_b, col_type and col_abbrev mentioned in the inputs description above.
collate_interactions_from_both_direction <- function (input_table, gene_a, gene_b, 
													  interaction_type_forward, interaction_type_reverse, 
													  interaction_abbrev_forward, interaction_abberv_reverse,
													  col_a, col_b, col_type, col_abbrev) {
	num_rows <- length(input_table[,1])
	
	new_table <- as.data.frame( matrix(ncol=4, nrow=2*num_rows))
	
	colnames( new_table) <- c(col_a, col_b, col_type, col_abbrev)
	
	new_table[1:num_rows, col_a] <- as.character(input_table[, gene_a] )
	new_table[1:num_rows, col_b] <- as.character(input_table[, gene_b] )
	new_table[1:num_rows, col_type] <- rep(  interaction_type_forward, num_rows)
	new_table[1:num_rows, col_abbrev] <- rep(  interaction_abbrev_forward, num_rows)
	
	new_table[(num_rows+1):(2*num_rows), col_a] <- as.character(input_table[, gene_b] )
	new_table[(num_rows+1):(2*num_rows), col_b] <- as.character(input_table[, gene_a] )
	new_table[(num_rows+1):(2*num_rows), col_type] <- rep( interaction_type_reverse, num_rows)
	new_table[(num_rows+1):(2*num_rows), col_abbrev] <- rep( interaction_abberv_reverse, num_rows)	
	
	new_table <- unique(new_table)
	
	new_table <- new_table[order(new_table[,col_a], new_table[,col_b]),]
	
	return(new_table)
}


#########################################################
## Insert directed network into both directions into the edge list table, so the direction that table joins are done won't affect analyses
## The table must have 4 columns only, which are: gene_a, gene_b, interaction type, interaction type abbreviation
### Same as above but for tbl_df objects
## Inputs:
## input_table - must have 4 columns only, which are: gene_a, gene_b, interaction type, interaction type abbreviation
## gene_a - column name for input table that cotains gene A of the interaction
## gene_b - column name for input table that cotains gene B of the interaction
## interaction_type_forward - information on the type of interaction for A ---> B 
## interaction_type_reverse - information on the type of interaction for B<---> A
## interaction_abbrev_forward - abbereviated information on the type of interaction for A ---> B 
## interaction_abberv_reverse - abbereviated information on the type of interaction for B ---> A 
## col_a - column name for gene A in the output table
## col_b - column name for gene B in the output table
## col_type - column name for the interaction type in the output table
## col_abbrev - column name for the interaction type abbreviation in the output table


## Output table:
## Output table has 4 columns, denoted by col_a, col_b, col_type and col_abbrev mentioned in the inputs description above.
collate_interactions_from_both_direction_tbl_df <- function (input_table, gene_a, gene_b, 
															 interaction_type_forward, interaction_type_reverse, 
															 interaction_abbrev_forward, interaction_abbrev_b,
															 col_a, col_b, col_type, col_abbrev) {
	
	num_rows <- count(input_table)[[1]]
	
	input_table <- as.data.frame(input_table)
	
	#print ( colnames(input_table))
	temp1 <- data.frame(rep(  interaction_type_forward, num_rows))
	temp2 <- data.frame(rep(  interaction_abbrev_forward, num_rows))
	
	# print(head(input_table))
	# print ( colnames(input_table) )
	# print(gene_a)
	# print(gene_b)

	left_table  <- dplyr::tbl_df( cbind(input_table[, gene_a], 
								        input_table[, gene_b] , 
								        data.frame(rep(  interaction_type_forward, num_rows)), 
								        data.frame(rep(  interaction_abbrev_forward, num_rows))) ) 
	colnames(left_table ) <- c( col_a, col_b, col_type, col_abbrev)			
	
	right_table <- dplyr::tbl_df( cbind(input_table[, gene_b], 
									    input_table[, gene_a] , 
									    data.frame(rep(  interaction_type_reverse, num_rows)), 
									    data.frame(rep(  interaction_abbrev_b, num_rows))) ) 
	
	colnames(right_table ) <- c( col_a, col_b, col_type, col_abbrev)			
	
	
	new_table <- dplyr::bind_rows( left_table, right_table) %>% 
		distinct () %>% 
		arrange_( .dots=lapply (c(col_a, col_b), as.name))

	
	return(new_table)
}



#########################################################
# Union tables together (i.e. similar to doing lots of rbind, but since the table is dynamically created so it is more efficient)
# Input: List of tables which have the same number of columns
# Output: Table with all input table unioned together.
efficient_union_of_tables <- function( list_of_tables,  stringsAsFactors=FALSE) {
	
	# Check tables has the same number of columns
	# Count the number of rows in all the tables
	prev_num_of_columns <- 0
	curr_num_of_columns <- 0 
	total_num_of_rows <- 0 
	
	for( i in 1:length( list_of_tables) ) {
		
		curr_num_of_columns <-	ncol(list_of_tables[[i]])
		
		
		if ( i > 1 ) {
			
			if ( curr_num_of_columns != prev_num_of_columns) {
				error_message <- paste( "efficient_union_of_tables: Input table in list location", i ,"has", curr_num_of_columns, 
										"rows. This is not the same previous tables which has", prev_num_of_columns, "rows." )
				
				stop(error_message )
			}
			
		}
		
		# Save the number of columns
		prev_num_of_columns <- curr_num_of_columns
		
		# Tally the total number of rows
		total_num_of_rows <- total_num_of_rows + nrow(list_of_tables[[i]])
	}
	
	
	# Create a table containing the union of all the tables
	total_num_of_columns <- prev_num_of_columns
	
	new_table <- as.data.frame( matrix(ncol=total_num_of_columns, nrow=total_num_of_rows), stringsAsFactors=stringsAsFactors)
	
	colnames(new_table) <- colnames(list_of_tables[[i]])
	
	up_to_row <- 0
	
	for( i in 1:length( list_of_tables) ) {
		
		num_rows <- nrow(list_of_tables[[i]])
		
		new_table[(up_to_row+1):(up_to_row+num_rows),  ] <- list_of_tables[[i]]
		
		up_to_row <- up_to_row + num_rows
	}
	
	return ( new_table)
}



#########################################################
## Form the triplet motifs and remove any duplicates, regardless of the orientation of gene A and gene B
## Need to remember that if the column names clashes between the tables, you'll need to append '.x' or '.y' to the column names 
## in the input join statements
## Inputs: 
# table_a - the genetic interaction network table
# table_b - list of all other types of interactions (i.e. kinase-substrate, protein-protein, transcription factor - target gene )
# table_c - list of all other types of interactions (i.e. kinase-substrate, protein-protein, transcription factor - target gene )
# join_ac - the join conditions between table_a and table_c
# join_bc - the join conditions between table_b and table_c. Need to include extra suffix '.x' or '.y' create internally by the dply::inner_join function
# selected_columns - The columns selected from the combined table, after the joins are made. 
#                    Need to include extra '.x' or '.y' create internally by the dply::inner_join function.
# column_names - After the columns are selected, rename them into the names you choose.

## Example inputs:
# my_join_ac <- c( "query_oln_id_edited" = "oln_id_a")
# my_join_bc <- c( "array_oln_id_edited"= "oln_id_a", "oln_id_b" = "oln_id_b" )
# my_selected_columns <- c("query_oln_id_edited", "array_oln_id_edited", 
# 						 "oln_id_b",  "interaction_type_abbrev.x",
# 						 "interaction_type_abbrev.y", "genetic_interaction_score", "p_value", "std_dev")
# my_column_names <- c( "oln_id_a", 
# 					  "oln_id_b",
# 					  "oln_id_c",
# 					  "type_ac",
# 					  "type_bc",
# 					  "genetic_interaction_score",
# 					  "p_value",
# 					  "std_dev" )  

form_triplet_motifs <- function (table_a, table_b, table_c, join_ac, join_bc, selected_columns, column_names ) { 
	
  my_debug <- 0
  
	old_column_names <- colnames(table_a)
	
	table_a <- tbl_df(table_a)
	table_b <- tbl_df(table_b)
	table_c <- tbl_df(table_c)
	
	temp_merged_table <- dplyr::inner_join( table_a, table_b, by= join_ac )  
	
	if( my_debug == 1) { print ( paste( "count(temp_merged_table) =", count(temp_merged_table)[[1]]) ) }
	
	
	# print ( paste ( colnames(temp_merged_table) )) 
	
	temp_merged_table <- inner_join( temp_merged_table, table_c, by= join_bc )
	
	if( my_debug == 1) { print ( paste( "count(temp_merged_table) =", count(temp_merged_table)[[1]]) ) }
	
	
	# print ( paste ( colnames(temp_merged_table) )) 
	
	
	triplet_motifs<- dplyr::select( temp_merged_table, one_of(selected_columns ) )  %>%
			 dplyr::rename_( .dots=setNames(as.list(selected_columns), column_names ) )
	
	
	if( my_debug == 1) { print ( paste( "count(triplet_motifs) =", count(triplet_motifs)[[1]]) ) }
	
	
	a_gt_b <- filter_( triplet_motifs, paste(column_names[1], ">", column_names[2] )   )
	
	if( my_debug == 1) { print ( paste( "count(a_gt_b) =", count(a_gt_b)[[1]]) ) }
	
	updated_column_names <- c(column_names[2], column_names[1], column_names[3], column_names[5],  column_names[4])
	
	## I would like to switch to using standard evalualtion of names, because this allows more flexibility in table names
	
	b_gt_a <- filter_( triplet_motifs, paste(column_names[2], ">", column_names[1] )   ) %>%
		rename_( .dots=setNames(as.list(column_names[1:5]), updated_column_names ) ) %>%
		dplyr::select( one_of(column_names))

	if( my_debug == 1) { print ( paste( "count(b_gt_a) =",count(b_gt_a)[[1]]) )}
	
	
	triplet_motifs <- dplyr::union(a_gt_b, b_gt_a )  %>%
						dplyr::group_by_( .dots=lapply ( column_names[1:5], as.name) ) %>%
						dplyr::summarise_(max_gi_score =  interp(~max(var), var = as.name(column_names[6])) )  %>%
						dplyr::rename_(.dots=setNames(list("max_gi_score"),c(column_names[6])))  %>%
						dplyr::filter_(paste( column_names[1], "!=", column_names[3] )  ) %>%
						dplyr::filter_(paste( column_names[2], "!=", column_names[3] )  ) %>%
						dplyr::arrange_( .dots=lapply (column_names[1:3], as.name))   %>%
					 	dplyr::ungroup()
	

	return(triplet_motifs)
}


#########################################################
## Count the number of each type of triplet motif.
## The triplet motif cosists of nodes a, b, and c. 
## Interaction A-B is fixed as a genetic interaction.
## We are interested in the type of interactions for edge A-C and B-C (i.e. type_ac and type_bc).
## User can chose the name of the columns for 'type_ac' and 'type_bc'
## Input:
## triplet_motifs - Table in the following format
# Source: local data frame [6 x 6]
# Groups: oln_id_a, oln_id_b, oln_id_c, type_ac [4]
# 
#	  oln_id_a oln_id_b oln_id_c type_ac type_bc genetic_interaction_score
#	  (chr)    (chr)    (chr)   (chr)   (chr)                     (dbl)
# 1   YBR060C  YAL019W  YBR160W      ku      ku                   -0.2641
# 2   YBR060C  YAL019W  YBR160W      ku       p                   -0.2641
# 3   YBR060C  YAL019W  YBR160W       p      ku                   -0.2641
# 4   YBR060C  YAL019W  YBR160W       p       p                   -0.2641
# 5 YBR111W-A  YAR002W  YBR010W       p       p                   -0.2409
# 6 YBR111W-A  YAR002W  YPL169C       p       p                   -0.2409

count_triplet_motifs_custom <- function (triplet_motifs, type_ac, type_bc, total_count ) {
	
	a_ge_b <- filter_( triplet_motifs, paste(type_ac, ">=", type_bc) ) %>%
				group_by_( .dots=lapply ( c(type_ac, type_bc), as.name)  ) %>%
				summarise(counts=n())
	
	columns_to_select <- c(type_ac, type_bc, "counts")
	
	b_gt_a <- filter_( triplet_motifs, paste(type_bc, ">", type_ac) ) %>%
				group_by_( .dots=lapply ( c(type_ac, type_bc), as.name) ) %>%
				summarise(counts=n()) %>%
				dplyr::rename_( .dots=setNames(list(type_ac, type_bc), c(type_bc, type_ac) ) ) %>%
				select( one_of( columns_to_select ))

	triplet_motif_counts <- dplyr::union( a_ge_b, b_gt_a)  %>%
		group_by_( .dots=lapply ( c(type_ac, type_bc), as.name) ) %>%	
		summarise( temp_counts=sum(counts))  %>%
		dplyr::rename_( .dots=setNames(list("temp_counts"), c(total_count) ) ) %>%
		arrange_( .dots=c(type_ac, type_bc) )
	

	return( triplet_motif_counts)
}


#########################################################
## Count the number of each type of triplet motif.
## Assume that there are the columns type_ac and type_bc.
## Input:
## triplet_motifs - Table in the following format
# Source: local data frame [6 x 6]
# Groups: oln_id_a, oln_id_b, oln_id_c, type_ac [4]
# 
#	  oln_id_a oln_id_b oln_id_c type_ac type_bc genetic_interaction_score
#	  (chr)    (chr)    (chr)   (chr)   (chr)                     (dbl)
# 1   YBR060C  YAL019W  YBR160W      ku      ku                   -0.2641
# 2   YBR060C  YAL019W  YBR160W      ku       p                   -0.2641
# 3   YBR060C  YAL019W  YBR160W       p      ku                   -0.2641
# 4   YBR060C  YAL019W  YBR160W       p       p                   -0.2641
# 5 YBR111W-A  YAR002W  YBR010W       p       p                   -0.2409
# 6 YBR111W-A  YAR002W  YPL169C       p       p                   -0.2409

count_triplet_motifs <- function (triplet_motifs ) {

	a_ge_b <- filter( triplet_motifs, type_ac >= type_bc ) %>%
		group_by( type_ac, type_bc) %>%
		summarise(counts=n())
	
	b_gt_a <- filter( triplet_motifs, type_bc > type_ac ) %>%
		group_by( type_ac, type_bc) %>%
		summarise(counts=n()) %>%
		rename ( type_ac= type_bc, type_bc = type_ac) %>%
		select( one_of( c("type_ac", "type_bc", "counts") ))
	
	triplet_motif_counts <- dplyr::union( a_ge_b, b_gt_a)  %>%
		group_by( type_ac, type_bc) %>% 
		summarise(total_count=sum(counts)) %>%
		arrange( type_ac, type_bc)
	
	return( triplet_motif_counts)
	
}

#########################################################

# Transpose the outputput from count_triplet_motifs_custom
# Input:
# triplet_motif_counts - output from count_triplet_motifs_custom
# type_ac, type_bc, total_count, motif_type - column names from triplet_motif_counts

# Output: triplet_motif_counts
# Each row contains the number of each triplet motif. Each column is one type of triplet motif.

transpose_count_triplet_motifs_results <- function (triplet_motif_counts,  type_ac, type_bc, total_count, motif_type) {

	## Mutate: this concatenate the columns type_ac and type_bc into one column called 'motif_type'
	## Remove columns type_ac and type_bc, and put motif_type column into the first column, the rest of the columns contain the counts for each triplet motif, then transpose the table

	mutate_dots <- ~paste(type_ac, type_bc, sep="")
	triplet_motif_counts <- mutate_(triplet_motif_counts, .dots=setNames(list(mutate_dots), c(motif_type) ) )
	triplet_motif_counts <- triplet_motif_counts[, c(motif_type, total_count)]
	
	colnames_to_use <- as.vector ( t( triplet_motif_counts[, motif_type] ))
	
	triplet_motif_counts <- t( triplet_motif_counts[,  total_count])
	
	colnames( triplet_motif_counts) <- colnames_to_use
	
	return( triplet_motif_counts)
}


#########################################################
#### Code to randomize a graph, keeping the degree distribution the same
### Uses the first two columns of the input table to create a graph
## Inputs:
## edge_list_table - Edge list table 
## col_a           - Name of the first column (A) from the input 'graph-table' to use
## col_b           - Name of the second column (B) from the input 'graph-table' to use
## directed        - Is the network a directed network, default to TRUE
## num_iterations  - The number of times the two edges are randomly chosen to shuffle them. The input value is the a multiple of the number of 'edges' in the network. 
##                 - If NULL the number of rounds of rewiring equal to the number of edges in the network
##                 - Additional information can be accessed here: https://lists.nongnu.org/archive/html/igraph-help/2014-04/msg00005.html
## min_num_neighbours - The minimum number of edges of similar joint degree distribution of the nodes to choose from. 
## max_num_tries 	  - Maximum number of tries to find an edge that does not result in a self-loop or multiple edges.
### Output:
### Edge list of rewired graph

rewire_graph <- function(edge_list_table, col_a, col_b, directed=TRUE, loops=TRUE, num_iterations=NULL, min_num_neighbours=10, max_num_tries=10, mode="dd" ) {

	graph_object <- graph.data.frame(edge_list_table[, c(col_a,col_b)], directed=directed, vertices=NULL)
	
	if ( is.null(num_iterations)  ) {
		num_iterations<- ecount(graph_object) 
	}  else if ( is.numeric( num_iterations ) ) {
		num_iterations <- ecount(graph_object) * num_iterations
	} else {
		stop ( paste ( "Input for number of iterations (", num_iterations, ") not invalid ", sep=""))
	}
	
	print (paste( "rewire_graph: num_iterations = ", num_iterations, sep="" ) )
	
	if ( mode== "dd" ) {

		print ("rewire_graph: Mode degree distribution")
	
		new_graph <- graph_object %>%
					 rewire(keeping_degseq(loops = loops, niter = num_iterations))
	
		graph_rewired_table <- dplyr::tbl_df(igraph::as_data_frame(new_graph, what="edges"))	
			
		graph_rewired_table <- dplyr::rename_(graph_rewired_table, .dots=setNames(list("from", "to"), c(col_a, col_b) )) %>% as.data.frame()
		

	} else if ( mode == "jdd") {
		
		print ("rewire_graph: Mode joint degree distribution")
		
		
		graph_rewired_table <- 	rewire_jdd_graph(edge_list_table, col_a, col_b, 
												min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, 
												num_iterations=num_iterations,  directed=directed ) %>% as.data.frame()
		
	}
	
	return(graph_rewired_table)
}

#########################################################
### Function: clean_genetic_interactions_table
## Description: Clean the genetic interactions edge list. The value of column A is greater than column B, no repetitive edges.

## For duplicated / multiple edges, the maximum negative genetic interaction score is kept. 
## Inputs:
## graph_table - edge list table 
## col_a - name of the first column (A) from the input 'graph-table' to use
## col_b - name of the second column (B) from the input 'graph-table' to use
## score - name of the score column from the input 'graph-table' to use

## Output: 
# regularized_table - A table with the following columns:
# col_a
# col_b
# score 
clean_genetic_interactions_table <- function (graph_table,  col_a, col_b, score ) {

	graph_table <- dplyr::tbl_df(graph_table)
	
	columns_selected <- c(col_a, col_b, score)
		
	a_ge_b <- filter_(graph_table, paste( col_a, ">=", col_b) ) %>%
		dplyr::select( one_of(columns_selected)) 
	
	b_gt_a <- filter_(graph_table, paste( col_b, ">", col_a) ) 
		#print ( colnames( a_ge_b) )
		#print( colnames( b_gt_a ) )
		#print ( as.data.frame( b_gt_a[1,]) )
		print ( summarize(a_ge_b,   count = n())   )
	    print ( summarize(b_gt_a,   count = n())   )

	b_gt_a <-	rename_(b_gt_a, .dots=setNames(list(col_a, col_b), c(col_b, col_a) ) )  %>%
		dplyr::select( one_of(columns_selected))
	
		#print( colnames( b_gt_a ) )
		#print ( as.data.frame( b_gt_a[1,]) )
	
	summarise_dots <- interp(~max(var), var = as.name(score))  #  ~mean(counts) 
	
	regularized_table <- dplyr::union( a_ge_b, b_gt_a ) # %>% 
	regularized_table <- dplyr::group_by_(regularized_table, .dots=lapply ( c(col_a, col_b), as.name) ) # %>%
	regularized_table <- dplyr::summarise_(regularized_table, .dots=setNames(list(summarise_dots), c(score))  ) # %>%
	regularized_table <- dplyr::ungroup(regularized_table)
	
	#print ( summarise_dots)
	#print ( paste (colnames(regularized_table)))
	
	regularized_table <- dplyr::select(regularized_table, one_of(columns_selected))
	
	return ( regularized_table)

}

### Rewire edge list table associated with a score attribute
### The scores are only here as a place-holder, they don't actually sync with the randomization so cannot actually be used.
### In the version of igraph used, the edge score / attributes are corrupted upon randomization.
## Inputs:
## graph_table - edge list table 
## col_a - name of the first column (A) from the input 'graph-table' to use
## col_b - name of the second column (B) from the input 'graph-table' to use
## score - name of the score column from the input 'graph-table' to use
## directed - is the network a directed network, default to FALSE
## num_iterations - number of rounds to perform the rewiring, if NULL the number of rounds of rewiring equal to the number of nodes in the network times 10
## min_num_neighbours - The minimum number of edges of similar joint degree distribution of the nodes to choose from. 
## max_num_tries 	  - Maximum number of tries to find an edge that does not result in a self-loop or multiple edges.

## Output: 
# randomized_table_with_score - A table with the following columns:
# col_a
# col_b
# score 
rewire_genetic_interactions_table <- function (graph_table,  col_a, col_b, score, directed=FALSE, num_iterations=NULL,
											    min_num_neighbours=10, max_num_tries=10, mode="dd") {
	
	regularized_table <- clean_genetic_interactions_table(graph_table,  col_a, col_b, score ) 

	print( paste( "rewire_genetic_interactions_network = ", count(regularized_table ) , sep="" ) )  
	
	system.time ( graph_rewired_table <- rewire_graph(regularized_table, col_a, col_b, directed=directed, loops=FALSE, num_iterations=num_iterations,
													  min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, mode=mode) )

	 randomized_table_with_score <- bind_cols(graph_rewired_table, regularized_table[,score] )
	
	# randomized_table_with_score <- graph_rewired_table
	
	return(randomized_table_with_score )
}

#########################################################

## Function: clean_graph_table
# Description: Clean the edge list table, where the value of column A is greater than column B. Remove duplicated edges.
## Inputs:
# graph_table - edge list table 
# col_a - name of the first column (A) from the input 'graph-table' to use
# col_b - name of the second column (B) from the input 'graph-table' to use
# directed - is the network a directed network, default to TRUE
## Output:
# regularized_table - A table with the following columns:
# col_a, col_b
clean_graph_table <- function (graph_table,  col_a, col_b, directed=FALSE ) {

	graph_table <- dplyr::tbl_df(graph_table)
	regularized_table <- graph_table
	
	if ( directed==FALSE) { 
		a_ge_b <- filter_(graph_table, paste( col_a, ">=", col_b) ) %>% 
			dplyr::select( one_of( c(col_a, col_b) ))
			
		b_gt_a <- filter_(graph_table, paste( col_b, ">", col_a) ) %>% 
			rename_( .dots=setNames(list(col_a, col_b), c(col_b, col_a) ) )  %>% 
			dplyr::select( one_of( c(col_a, col_b) ))
		
		regularized_table <- dplyr::union( a_ge_b, b_gt_a ) # %>% 
		regularized_table <- dplyr::select(regularized_table, one_of( c(col_a, col_b) )) # %>%
		regularized_table <- dplyr::distinct(regularized_table ) 
	} else {
	
		regularized_table <- dplyr::distinct(graph_table ) 

	}
		
	return(regularized_table )

}


### Rewire graph by inputing an edge list table  
## Inputs:
## graph_table - edge list table 
## col_a - name of the first column (A) from the input 'graph-table' to use
## col_b - name of the second column (B) from the input 'graph-table' to use
## directed - is the network a directed network, default to TRUE
## num_iterations - number of rounds to perform the rewiring, if NULL the number of rounds of rewiring equal to the number of nodes in the network times 10

## Outputs:
# graph_rewired_table - A table with the following columns:
# col_a, col_b
rewire_graph_table <- function (graph_table,  col_a, col_b, directed=TRUE, loops=TRUE, num_iterations=NULL, 
								min_num_neighbours=10, max_num_tries=10, mode="dd") {
	
	regularized_table <- clean_graph_table (graph_table,  col_a, col_b, directed=directed ) 

	system.time ( graph_rewired_table <- rewire_graph(regularized_table, col_a, col_b, directed=directed, loops=loops, num_iterations=num_iterations,
													  min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, mode=mode) )

	return(graph_rewired_table )
}


#########################################################
### Get the kinase-substrate network, transcription factor-target network, protein-protein interaction network
### Shuffle each one and then combine them together into one table.

## Inputs: 
## Requires the following input tables, with these specific columns of the same name.
## 		* kinase_network must have these columns: kinase_oln_id, target_oln_id
##  	* sbi_interactome must have these columns: oln_id_a, oln_id_b
##  	* rewired_tf_network must have these columns: regulator_oln_id, target_oln_id
##  num_iteractions: number of times to rewire the network. A higher number of times ensure that the network is well randomized. Default 1000 times

## Outputs:
## Table of all interactions, with the following columns
## oln_id_a - ordered locus name (OLN) of gene A
## oln_id_b - OLN of gene B
## interaction_type - description of the type of interactions and the direction of interactions. The kinase-substrate network and transcription factor-target network have up or down direction. Protein-protein interaction network do not have have direction.
## interaction_type_abbrev - kinase-up (ku), kinase-down (kd), transcription factor-up (tu), transcription factor-down (td), protein-protein interaction (p)


rewired_interaction_network <- function (kinase_network, sbi_interactome, tf_network, num_iterations = NULL,  min_num_neighbours=10, max_num_tries=10, mode="dd" ) {
	
	rewired_kinase_network <-  shuffle_table(kinase_network) 
	
	print( paste( "rewired_kinase_network = ", count(rewired_kinase_network ) , sep="" ) )  

    # print(colnames(rewired_kinase_network ))
	
	rewired_kinase_network <- rewire_graph_table( rewired_kinase_network, "kinase_oln_id", "target_oln_id", directed=TRUE, num_iterations=num_iterations,
							min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, mode=mode ) 
	
	# print( paste( "rewired_kinase_network = ", count(rewired_kinase_network ) , sep="" ) )  

	rewired_kinase_network <-	collate_interactions_from_both_direction_tbl_df(  rewired_kinase_network, "kinase_oln_id", "target_oln_id", 
														  "kinase-substrate down", "kinase-substrate up",
														  "kd", "ku", 
														  "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev") 
	
	# print( paste( "rewired_kinase_network = ", count(rewired_kinase_network ) , sep="" ) )  
	
	rewired_ppi_network <- 	shuffle_table(sbi_interactome) 
	
    print( paste( "rewired_ppi_network = ", count(rewired_ppi_network ) , sep="" ) )  

	
	rewired_ppi_network <-	rewire_graph_table( rewired_ppi_network,  "oln_id_a", "oln_id_b", directed=FALSE, 
							num_iterations=num_iterations,
							 min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, mode=mode ) 
							 
	# print( paste( "rewired_ppi_network = ", count(rewired_ppi_network ) , sep="" ) )  

	rewired_ppi_network <-	collate_interactions_from_both_direction_tbl_df( rewired_ppi_network, "oln_id_a", "oln_id_b", 
														  "protein-protein", "protein-protein",
														  "p", "p", 
														  "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev") 
	
	# print( paste( "rewired_ppi_network = ", count(rewired_ppi_network ) , sep="" ) )  

	
	rewired_tf_network <- shuffle_table(tf_network) 
	
	print( paste( "rewired_tf_network = ", count(rewired_tf_network ) , sep="" ) )  

	
	rewired_tf_network <-	rewire_graph_table( rewired_tf_network, "regulator_oln_id", "target_oln_id", directed=TRUE, num_iterations=num_iterations,
							 min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, mode=mode )  
		
	# print( paste( "rewired_tf_network = ", count(rewired_tf_network ) , sep="" ) )  

		
	rewired_tf_network <-	collate_interactions_from_both_direction_tbl_df( rewired_tf_network, "regulator_oln_id", "target_oln_id", 
														 "transcription factor-target down", "transcription factor-target up",
														 "td", "tu", 
														 "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev")
	
	# print( paste( "rewired_tf_network = ", count(rewired_tf_network ) , sep="" ) )  

	
	rewired_interactions_combined <- as.data.frame( dplyr::bind_rows( rewired_tf_network, rewired_kinase_network, rewired_ppi_network) )
	
	return( rewired_interactions_combined)
	
}


### Get the kinase-substrate network, transcription factor-target network, protein-protein interaction network
### Shuffle each network.
### Get the driver nodes separately for the kinase-substrate network and the transcription factor-target network.
### Combine all the networks together into one table.

## Inputs: 
## Requires the following input tables, with these specific columns of the same name.
## 		* kinase_network must have these columns: kinase_oln_id, target_oln_id
##  	* sbi_interactome must have these columns: oln_id_a, oln_id_b
##  	* rewired_tf_network must have these columns: regulator_oln_id, target_oln_id
##  num_iteractions: number of times to rewire the network. A higher number of times ensure that the network is well randomized. Default 1000 times
## remove.multiple   = FALSE, - remove duplicated edges  
## remove.loops      = TRUE,  - remove self-loops
## delete_leaf_nodes = TRUE,  - If TRUE, do not consider nodes which have not out degrees in the controllability calculations. If FALSE, consider all nodes.

## Outputs:
## Table of all interactions, with the following columns
## oln_id_a - ordered locus name (OLN) of gene A
## oln_id_b - OLN of gene B
## interaction_type - description of the type of interactions and the direction of interactions. The kinase-substrate network and transcription factor-target network have up or down direction. Protein-protein interaction network do not have have direction.
## interaction_type_abbrev - kinase-up (ku), kinase-down (kd), transcription factor-up (tu), transcription factor-down (td), protein-protein interaction (p)

controllability_of_rewired_network <- function (kinase_network, sbi_interactome, tf_network, num_iterations = NULL, 
												remove.multiple = FALSE, remove.loops=TRUE, delete_leaf_nodes=TRUE ) {
	
	rewired_kinase_network <-  shuffle_table(kinase_network) %>% 
		rewire_graph_table( "kinase_oln_id", "target_oln_id", directed=TRUE, num_iterations=num_iterations ) 

	rewired_kinase_network_collated <-	collate_interactions_from_both_direction_tbl_df( rewired_kinase_network, "kinase_oln_id", "target_oln_id", 
														  "kinase-substrate down", "kinase-substrate up",
														  "kd", "ku", 
														  "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev") 
	
	rewired_ppi_network <- 	shuffle_table(sbi_interactome) %>% 
		rewire_graph_table(  "oln_id_a", "oln_id_b", directed=FALSE, num_iterations=num_iterations ) %>%
		collate_interactions_from_both_direction_tbl_df(  "oln_id_a", "oln_id_b", 
														  "protein-protein", "protein-protein",
														  "p", "p", 
														  "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev") 
	
	rewired_tf_network <- shuffle_table(tf_network) %>%  
		rewire_graph_table(  "regulator_oln_id", "target_oln_id", directed=TRUE, num_iterations=num_iterations ) 
	
	rewired_tf_network_collated <- collate_interactions_from_both_direction_tbl_df( rewired_tf_network, "regulator_oln_id", "target_oln_id", 
														 "transcription factor-target down", "transcription factor-target up",
														 "td", "tu", 
														 "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev")
	
	rewired_interactions_combined <- dplyr::bind_rows( rewired_tf_network_collated, rewired_kinase_network_collated, rewired_ppi_network)
	
	
	# Controllability analysis
	
	kinase_network_driver_nodes <- controllability_analysis ( rewired_kinase_network, "kinase_oln_id", "target_oln_id", 
															  remove.multiple =remove.multiple, remove.loops=remove.loops, 
															  delete_leaf_nodes=delete_leaf_nodes  )  
	
	tf_network_driver_nodes     <- controllability_analysis ( rewired_tf_network, "regulator_oln_id", "target_oln_id", 
															  remove.multiple =remove.multiple, remove.loops=remove.loops, 
															  delete_leaf_nodes=delete_leaf_nodes  ) 
	
	
	return( list( rewired_network= rewired_interactions_combined, 
				  kinase_driver_nodes=kinase_network_driver_nodes,  
				  tf_driver_nodes=tf_network_driver_nodes) ) 
	
}




#########################################################
## Function: shuffle_table
## Shuffle the rows of the input table, without replacement.
## Input: A data frame table
## Output: 
## Returns the shuffled table. Cast to a 'tibble' data frame using the dplyr::tbl_df function. 
## Please refer to the tibble and dplyr libraries for more information.

shuffle_table <- function (input_table) {
	
	input_table <- tbl_df(input_table)
	
	shuffled_indicies <- sample(1:count(input_table)[[1]], count(input_table)[[1]], replace = FALSE)
	
	shuffled_table <- input_table[shuffled_indicies,]
	
	return(shuffled_table)
}

#########################################################

## This is to skip a few random number so that it will be different next time things are run using mclapply
	skip.streams <- function(n) {
		x <- .Random.seed
		for (i in seq_len(n))
			x <- nextRNGStream(x)
		assign('.Random.seed', x, pos=.GlobalEnv)
	}


#########################################################

### Take a list of results from applying multiple rounds of 'run_one_randomized_trial'
# and merge them into one table (i.e. the reduce step of map-reduce)
# The motif type is in the first column, the rest of the columns contain the counts for each randomized trial.
concat_motif_counts_list_into_table <- function (list_of_randomized_triplet_motif_counts ) {
	
	temp_table <- as.data.frame( list_of_randomized_triplet_motif_counts[[1]] )
	vector_column_names <- c("type_ac", "type_bc", "1" )
	colnames( temp_table) <- vector_column_names
	
	if (length(list_of_randomized_triplet_motif_counts) > 1 ) { 
		# Start renaming from second element in the list 
		for ( i in 2:length(list_of_randomized_triplet_motif_counts)) {
	
				temp_table <- dplyr::full_join ( temp_table, as.data.frame(list_of_randomized_triplet_motif_counts[[i]]), by=c("type_ac", "type_bc"))
				
				## Need to keep renaming the columns 2 to i, otherwise the join will not work properly as it gets the 'my_total_count.x' and 'my_total_count.y'
				## from the table joins confused. 
				vector_column_names  <- c(vector_column_names , as.character(i))
				colnames( temp_table) <- vector_column_names
		}
	}
		
	## Mutate: this concatenate the columns type_ac and type_bc into one column called 'motif_type'
	temp_table <- mutate(temp_table, motif_type=paste(type_ac, type_bc, sep=""))
	
	## Remove columns type_ac and type_bc, and put motif_type column into the first column, the rest of the columns contain the counts for each triplet motif
	# temp_table <- temp_table[, c("motif_type", as.character(1:length(list_of_randomized_triplet_motif_counts)))]

	
	colnames_to_use <- as.vector ( t( temp_table[, "motif_type"] ))
	
	temp_table <- t( temp_table[,  as.character(1:length(list_of_randomized_triplet_motif_counts))]) 
	
	colnames( temp_table) <- colnames_to_use
	
	
	return ( temp_table)
}



#########################################################

# Given the observed and randomization results, calculate the bootstrap p-values
# Inputs:
# observed_counts_table - table of observed counts, each triplet motif is represented by a column
# randomized_counts_table <- table of randomization counts, each triplet motif is represented by a column, counts for each randomized trial is represented by a row

# Output:
# A list containting the following elements
# $enriched - a vector of the p-values from testing for enrichment, names represent triplet motifs
# $depleted - a vector of the p-values from testing for depletion, names represent triplet motifs

calculate_boostrap_p_values <- function( observed_counts_table, randomized_counts_table ) {

	observed_counts_table <- tbl_df(as.data.frame( observed_counts_table) )
	randomized_counts_table <- tbl_df(as.data.frame ( randomized_counts_table) )
	num_randomization_trials <- nrow(randomized_counts_table)
	
	list_of_columns_to_test <- unique ( c(colnames( observed_counts_table), colnames(randomized_counts_table) ) )
	
	obs_gt_random <- rep(NA, length(list_of_columns_to_test))
	obs_lt_random <- rep(NA, length(list_of_columns_to_test))
	
	names( obs_gt_random) <- list_of_columns_to_test
	names( obs_lt_random) <- list_of_columns_to_test
	
	
	for  ( i in list_of_columns_to_test ) {
		
		# Deal with case where the triplet motif is not observed in the real data
		observed_value <- 0 
		if ( is.element(i,  colnames( observed_counts_table) ) ) {
			
			observed_value <- observed_counts_table[[1,i]]
		}
		
		if ( is.element(i,  colnames( randomized_counts_table) ) ) {
			obs_gt_random[i] <- length( which( observed_value > randomized_counts_table[,i]) ) 
			obs_lt_random[i] <- length( which( observed_value < randomized_counts_table[,i]) ) 
		} else {
			obs_gt_random[i] <- num_randomization_trials
			obs_lt_random[i] <- 0
		}	
			
	}
	
	obs_gt_random <- 1 - obs_gt_random/num_randomization_trials
	obs_lt_random <- 1 - obs_lt_random/num_randomization_trials
	
	return ( list ( enriched=obs_gt_random, depleted=obs_lt_random ))
}

#########################################################

## Get the randomized p-values, bonferroni adjusted p-values, median, mean and sd of the triplet motifs counts of all the randomized networks
## Test for enrichment as well as depletion. Use two-sided p-value
## Inputs: 
##  observed_counts_table - table of observed counts
##  randomized_counts_table - table of counts from many randomized networks
##
## Outputs:
## Table containing the following columns:
##    motif_type, observed_counts, 
##    enrichment_raw_p_values, enrichment_adj_p_values, 
##    depletion_raw_p_values, depletion_adj_p_values, 
##    mean, sd, median, ratio


get_full_results_table <- function (observed_counts_table, randomized_counts_table, p.value=0.05) {
	
	randomized_counts_table[is.na(randomized_counts_table)] <- 0 
	
	observed_counts_table <- tbl_df(as.data.frame( observed_counts_table) )
	randomized_counts_table <- tbl_df(as.data.frame ( randomized_counts_table) )
	
	## Calculate bootstrap p-value
	bootstrap_p_values <- calculate_boostrap_p_values( observed_counts_table, randomized_counts_table ) 
	
	## Perform bonferroni corrections
	bootstrap_p_values_enriched_updated <-  sapply ( (bootstrap_p_values$enriched)*length(bootstrap_p_values$enriched), 
													 function(x) ifelse ( x < p.value/2, x, NA)  )
	bootstrap_p_values_enriched_updated <- merge ( bootstrap_p_values$enriched, bootstrap_p_values_enriched_updated, by="row.names")
	colnames(bootstrap_p_values_enriched_updated ) <- c("motif_type", "enrichment_raw_p_values", "enrichment_adj_p_values")
	
	
	bootstrap_p_values_depleted_updated <-  sapply ( ( bootstrap_p_values$depleted)*length(bootstrap_p_values$depleted), 
													 function(x) ifelse ( x < p.value/2, x, NA)  )
	bootstrap_p_values_depleted_updated <- merge ( bootstrap_p_values$depleted, bootstrap_p_values_depleted_updated, by="row.names")
	colnames(bootstrap_p_values_depleted_updated ) <- c("motif_type", "depletion_raw_p_values", "depletion_adj_p_values")
	
	# Add motif_type column name to observed counts table
	observed_counts_table_updated <- as.data.frame(t(observed_counts_table)) 
	colnames( observed_counts_table_updated ) <- c("observed_counts")
	observed_counts_table_updated <- tibble::rownames_to_column(observed_counts_table_updated, "motif_type")
	
	# Calculate the mean and sd of triplet motifs counts among the randomized network
	bootstrap_mean <- as.data.frame(t( round(summarise_each ( randomized_counts_table, funs(mean) ) ) )) %>% tibble::rownames_to_column("motif_type")
	colnames( bootstrap_mean)[2] <- "mean"
	
	bootstrap_sd   <- as.data.frame(t(  round(summarise_each ( randomized_counts_table, funs(sd) ) ) ) ) %>% tibble::rownames_to_column("motif_type")
	colnames( bootstrap_sd)[2] <- "sd"
	
	bootstrap_median   <- as.data.frame(t(  round(summarise_each ( randomized_counts_table, funs(median) ) ) ) ) %>% tibble::rownames_to_column("motif_type")
	colnames( bootstrap_median)[2] <- "median"
	
	bootstrap_mad   <- as.data.frame(t(  round(summarise_each ( randomized_counts_table, funs(mad) ) ) ) ) %>% tibble::rownames_to_column("motif_type")
	colnames( bootstrap_mad)[2] <- "mad"
	
	# Count how many values are not NA
	my_count_values <- function(x) {  length(which(!is.na(x)))}
	
	bootstrap_count   <- as.data.frame(t(  summarise_each ( randomized_counts_table, funs(my_count_values)  ) ) )  %>% 
						 tibble::rownames_to_column("motif_type")
	colnames( bootstrap_count)[2] <- "count"

	# Join all the able together
	final_table_joined <- full_join ( observed_counts_table_updated, bootstrap_p_values_enriched_updated, by="motif_type") %>%
		full_join ( bootstrap_p_values_depleted_updated, by="motif_type")  %>%
		full_join ( bootstrap_mean, by="motif_type")  %>%
		full_join ( bootstrap_sd, by="motif_type")  %>%
		full_join ( bootstrap_median, by="motif_type")  %>%
		full_join ( bootstrap_mad, by="motif_type")  %>%
		full_join (bootstrap_count, by="motif_type")  %>%
		arrange(motif_type) 
	
	# Some triplet motifs do not have observed counts and only appear in full outer join
	final_table_joined[ is.na(final_table_joined[, "observed_counts"]), "observed_counts"]	 <- 0
	
	final_table_joined <-	mutate(final_table_joined, ratio=round(observed_counts/mean,2) )
	
	return( final_table_joined )
}


#########################################################

### Function: run_one_randomized_trial_compare_with_dataset
### Perform one round of randomization of all the networks. Form the triplet motifs and compare a set of query genes against the triplet motifs.
### Counts for protein A and B are combined together, as they are symmetrical in the triplet motif. Counts for protein C by itself.

## Get the randomized p-values, bonferroni adjusted p-values, median, mean and sd from comparing randomized triplet motifs dataset with a data table. 
## Tabulate results for all the randomized networks.
## Test for enrichment as well as depletion. Use two-sided p-value
## Inputs: 
## x: the variable x does nothing, just a place holder for the iterations in mclapply 
## Requires the following input tables, with these specific columns of the same name.
## 		* filtered_costanzo_stringent table must have these columns: oln_id_a, oln_id_b, genetic_interation_score
## 		* kinase_network must have these columns: kinase_oln_id, target_oln_id
##  	* sbi_interactome must have these columns: oln_id_a, oln_id_b
##  	* rewired_tf_network must have these columns: regulator_oln_id, target_oln_id

## FUNCT: Requires a function for comparing the triplet motifs with the query dataset. This function contain the rule for comparing the data tables.
## 		  An example function that can be used is 'count_triplet_motifs_against_gene_list'
## ...: Additional tables that are co-analyzed with the triplet motifs
## Outputs:
# A list of outputs in the same format as the output of the input function FUNCT.

run_one_randomized_trial_compare_with_dataset <- function (x, filtered_costanzo_stringent, kinase_network, sbi_interactome, tf_network, 
														   num_iterations=NULL, FUNCT,
														   min_num_neighbours=10, max_num_tries=10, mode="dd", use_rewired_interactions=FALSE,
														   ...) {
	
	rewired_filtered_costanzo_stringent <- rewire_genetic_interactions_table(filtered_costanzo_stringent, "oln_id_a",
																			 "oln_id_b", "genetic_interaction_score", 
																			 directed=FALSE, num_iterations=num_iterations,
																			 min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, mode=mode)
	
	rewired_interactions_combined <- rewired_interaction_network(kinase_network, sbi_interactome, tf_network, num_iterations=num_iterations,
																 min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, mode=mode)
	
	my_join_ac <- c( "oln_id_a" = "oln_id_a")
	my_join_bc <- c( "oln_id_b.x"= "oln_id_a", "oln_id_b.y" = "oln_id_b" )
	my_selected_columns <- c("oln_id_a", "oln_id_b.x", 
							 "oln_id_b.y",  "interaction_type_abbrev.x",
							 "interaction_type_abbrev.y", "genetic_interaction_score")
	my_column_names <- c( "oln_id_a", 
						  "oln_id_b",
						  "oln_id_c",
						  "type_ac",
						  "type_bc",
						  "genetic_interaction_score" )      
	
	rewired_triplet_motifs_costanzo <- form_triplet_motifs(rewired_filtered_costanzo_stringent, rewired_interactions_combined, 
														   rewired_interactions_combined,
														   my_join_ac, my_join_bc, my_selected_columns, my_column_names)
	
	## Counts of motifs in rewired network, filtering on the dataset provided as input as comparison. 
	
	if ( use_rewired_interactions == FALSE){
		rewired_counts_motifs_compare_with_dataset <- FUNCT ( rewired_triplet_motifs_costanzo, ...) 
	} else {
		rewired_counts_motifs_compare_with_dataset <- FUNCT ( rewired_triplet_motifs_costanzo, rewired_interactions_combined, ...) 
	}
	
	return( rewired_counts_motifs_compare_with_dataset )
	
}


## Function: run_one_randomized_trial_controllability
### Perform one round of randomization of all the networks. Form the triplet motifs. 
### Find the driver nodes in the randomized kinase-substrate network.
### Find the driver nodes in the randomized transcription factor-target gene network.
### Count the number of each type of triplet motifs that matches the kinase driver nodes.
### Count the number of each type of triplet motifs that matches the transcription factor driver nodes.
### Counts for protein A and B are combined together, as they are symmetrical in the triplet motif. Counts for protein C by itself.

## Get the randomized p-values, bonferroni adjusted p-values, median, mean and sd from comparing randomized triplet motifs dataset with a data table. 
## Tabulate results for all the randomized networks.
## Test for enrichment as well as depletion. Use two-sided p-value
## Inputs: 
## x: the variable x does nothing, just a place holder for the iterations in mclapply 
## Requires the following input tables, with these specific columns of the same name.
## 		* filtered_costanzo_stringent table must have these columns: oln_id_a, oln_id_b, genetic_interation_score
## 		* kinase_network must have these columns: kinase_oln_id, target_oln_id
##  	* sbi_interactome must have these columns: oln_id_a, oln_id_b
##  	* rewired_tf_network must have these columns: regulator_oln_id, target_oln_id

## FUNCT: Requires a function for comparing the triplet motifs with the query dataset. This function contain the rule for comparing the data tables.
## 		  An example function that can be used is 'count_triplet_motifs_against_gene_list'
## ...: Additional tables that are co-analyzed with the triplet motifs
## Outputs:
# A list of outputs in the same format as the output of the input function FUNCT.


run_one_randomized_trial_controllability <- function (x, filtered_costanzo_stringent, kinase_network, sbi_interactome, tf_network, 
														   num_iterations=NULL, 
													       remove.multiple = FALSE, remove.loops=TRUE, delete_leaf_nodes=TRUE, 
													       FUNCT, ...) {
	
	rewired_filtered_costanzo_stringent <- rewire_genetic_interactions_table(filtered_costanzo_stringent, "oln_id_a",
																			"oln_id_b", "genetic_interaction_score", 
																			directed=FALSE, num_iterations=num_iterations)
	
	rewired_network_results <- controllability_of_rewired_network(kinase_network, sbi_interactome, tf_network, num_iterations=num_iterations, 
																  remove.multiple = remove.multiple, remove.loops=remove.loops, 
																  delete_leaf_nodes=delete_leaf_nodes)
	
	kinase_driver_nodes <- rewired_network_results$kinase_driver_nodes
	tf_driver_nodes <- rewired_network_results$tf_driver_nodes
	
	rewired_interactions_combined <-rewired_network_results$rewired_network 
	
	my_join_ac <- c( "oln_id_a" = "oln_id_a")
	my_join_bc <- c( "oln_id_b.x"= "oln_id_a", "oln_id_b.y" = "oln_id_b" )
	my_selected_columns <- c("oln_id_a", "oln_id_b.x", 
							 "oln_id_b.y",  "interaction_type_abbrev.x",
							 "interaction_type_abbrev.y", "genetic_interaction_score")
	my_column_names <- c( "oln_id_a", 
						  "oln_id_b",
						  "oln_id_c",
						  "type_ac",
						  "type_bc",
						  "genetic_interaction_score" )      
	
	rewired_triplet_motifs_costanzo <- form_triplet_motifs(rewired_filtered_costanzo_stringent, rewired_interactions_combined, 
														   rewired_interactions_combined,
														   my_join_ac, my_join_bc, my_selected_columns, my_column_names)
	
	## Counts of motifs in rewired network, filtering on the dataset provided as input as comparison. 
	rewired_counts_motifs_compare_kinase_drivers <- FUNCT ( rewired_triplet_motifs_costanzo, kinase_driver_nodes, ...) 
	
	rewired_counts_motifs_compare_tf_drivers <- FUNCT ( rewired_triplet_motifs_costanzo, tf_driver_nodes, ...) 
	
	return( list ( kinase_driver_node_counts=rewired_counts_motifs_compare_kinase_drivers, 
				   tf_driver_node_counts=rewired_counts_motifs_compare_tf_drivers) )
	
}




run_one_randomized_trial_random_edges <- function (x, filtered_costanzo_stringent, kinase_network, sbi_interactome, tf_network, 
														   num_iterations=NULL, FUNCT,
														   min_num_neighbours=10, max_num_tries=10, mode="dd",
														   ...) {
	
	rewired_filtered_costanzo_stringent <- rewire_genetic_interactions_table(filtered_costanzo_stringent, "oln_id_a",
																			"oln_id_b", "genetic_interaction_score", 
																			directed=FALSE, num_iterations=num_iterations,
																			min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, mode=mode)
	
	rewired_interactions_combined <- rewired_interaction_network(kinase_network, sbi_interactome, tf_network, num_iterations=num_iterations,
																 min_num_neighbours=min_num_neighbours, max_num_tries=max_num_tries, mode=mode)
	
	
	my_join_ac <- c( "oln_id_a" = "oln_id_a")
	my_join_bc <- c( "oln_id_b.x"= "oln_id_a", "oln_id_b.y" = "oln_id_b" )
	my_selected_columns <- c("oln_id_a", "oln_id_b.x", 
							 "oln_id_b.y",  "interaction_type_abbrev.x",
							 "interaction_type_abbrev.y", "genetic_interaction_score")
	my_column_names <- c( "oln_id_a", 
						  "oln_id_b",
						  "oln_id_c",
						  "type_ac",
						  "type_bc",
						  "genetic_interaction_score" )      
	
	rewired_triplet_motifs_costanzo <- form_triplet_motifs(rewired_filtered_costanzo_stringent, rewired_interactions_combined, 
														   rewired_interactions_combined,
														   my_join_ac, my_join_bc, my_selected_columns, my_column_names)
	
	## Counts of motifs in rewired network, filtering on the dataset provided as input as comparison. 
	rewired_counts_motifs_compare_with_dataset <- FUNCT ( rewired_triplet_motifs_costanzo, ...) 
	
	return( rewired_counts_motifs_compare_with_dataset )
	
}



#########################################################

## Functipn: count_triplet_motifs_against_gene_list
## Description: Count the number of essential genes in each triplet motif 
## Inputs:
## triplet_motifs_list: A table containing the global list of triplet motifs
##	This input table must contain the following columns, 
#   oln_id_a - The Ordered Locus Name of Gene A
#   oln_id_b - The Ordered Locus Name of Gene B
#   oln_id_c - The Ordered Locus Name of Gene C
#   type_ac - The type of interactions between gene A and gene C, one of p, ku, kd, tu, td
#   type_bc - The type of interactions between gene B and gene C, one of p, ku, kd, tu, td  
## gene_list: A vector containig a list of essential genes in Ordered Locus Names.
###           This list could also denote any other type of genes e.g. mutants with a particular phenotype etc..)
## Outputs:
## A table containing the following columns:
## type_ac: Type of interaction between gene A and gene C of triplet motif
## type_bc: Type of interaction between gene B and gene C of triplet motif
## count_a_and_b: Total number of motifs with essential gene in either gene A or gene B of the triplet motif. 
## count_c: Total number of motifs with essential gene in gene C of the triplet motif, I will count this as 1
count_triplet_motifs_against_gene_list <- function(triplet_motifs_list, gene_list) {
	
	## Make sure triplet_motifs_list is ungrouped
	triplet_motifs_list <- ungroup(triplet_motifs_list)
	
	## Join the table for essential genes, so that we can count the number of essential genes in each column
	gene_list_edited <- as.data.frame( cbind( oln_id=gene_list, is_in_gene_list=rep(1, length(gene_list[,1]) ) ))
	
	triplet_motif_and_gene_lists <- left_join ( triplet_motifs_list, gene_list_edited, by=c("oln_id_a" = "oln_id"))  %>% 
		rename( oln_id_a_attribute=is_in_gene_list) %>%
		left_join (  gene_list_edited, by=c("oln_id_b" = "oln_id")) %>% 
		rename( oln_id_b_attribute=is_in_gene_list) %>%
		left_join (  gene_list_edited, by=c("oln_id_c" = "oln_id")) %>% 
		rename( oln_id_c_attribute=is_in_gene_list) 
	
	## count_a_and_b: If the attribute is found in either gene A or gene B of the triplet motif, I will count this as 1
	## count_c: If the attribute is found in gene C of the triplet motif, I will count this as 1
	

	### Do this for protein C
	a_ge_b <- filter( triplet_motif_and_gene_lists, type_ac >= type_bc ) %>%
		dplyr::distinct( type_ac, type_bc, oln_id_c, oln_id_c_attribute) %>%
		dplyr::group_by( type_ac, type_bc) %>%
		dplyr::summarise( count_c = sum( ifelse(is.na(oln_id_c_attribute),0,1)))
	
	b_gt_a <- filter( triplet_motif_and_gene_lists, type_bc > type_ac ) %>%
		dplyr::distinct( type_ac, type_bc, oln_id_c, oln_id_c_attribute) %>%
		dplyr::group_by( type_ac, type_bc) %>%
		dplyr::summarise(  count_c = sum( ifelse(is.na(oln_id_c_attribute),0,1)))  %>%
		dplyr::rename ( type_ac= type_bc, type_bc = type_ac) 
	
	count_genes_in_motifs_protein_c <- dplyr::union( a_ge_b, b_gt_a)  %>%
		group_by( type_ac, type_bc) %>% 
		summarise(count_c=sum(count_c)) %>%
		arrange( type_ac, type_bc)
	
	
	### Do this for protein A and B
	a_ge_b_AB <- filter( triplet_motif_and_gene_lists, type_ac >= type_bc ) %>%
		dplyr::distinct( type_ac, type_bc, oln_id_a, oln_id_b, oln_id_a_attribute, oln_id_b_attribute) %>%
		dplyr::group_by( type_ac, type_bc) %>%
		dplyr::summarise( count_a_and_b=sum(ifelse(is.na(oln_id_a_attribute) & is.na(oln_id_b_attribute),0,1)) )
	
	b_gt_a_AB <- filter( triplet_motif_and_gene_lists, type_bc > type_ac ) %>%
		dplyr::distinct( type_ac, type_bc, oln_id_a, oln_id_b, oln_id_a_attribute, oln_id_b_attribute) %>%
		dplyr::group_by( type_ac, type_bc) %>%
		dplyr::summarise( count_a_and_b=sum(ifelse(is.na(oln_id_a_attribute) & is.na(oln_id_b_attribute),0,1)) )  %>%
		dplyr::rename ( type_ac= type_bc, type_bc = type_ac) 
	
	count_genes_in_motifs_protein_a_and_b <- dplyr::union( a_ge_b_AB, b_gt_a_AB)  %>%
		group_by( type_ac, type_bc) %>% 
		summarise( count_a_and_b=sum(count_a_and_b) ) %>%
		arrange( type_ac, type_bc)
	
	
	count_genes_in_motifs <- full_join ( count_genes_in_motifs_protein_a_and_b, count_genes_in_motifs_protein_c, 
										 by=c("type_ac" = "type_ac", "type_bc" = "type_bc"))   %>%
							  arrange( type_ac, type_bc)
	
	return( count_genes_in_motifs)		
}


#########################################################


count_triplet_motifs_in_protein_complexes_helper <- function ( triplet_motifs_list, protein_complexes) {
	
	motif_and_protein_complexes <- left_join ( triplet_motifs_list, protein_complexes, by=c("oln_id_a" = "oln_id") ) %>%
		rename( complex_id_a = complex_id) %>%
		left_join ( protein_complexes, by=c("oln_id_b" = "oln_id") ) %>%
		rename( complex_id_b = complex_id) %>%
		left_join ( protein_complexes, by=c("oln_id_c" = "oln_id") ) %>%
		rename( complex_id_c = complex_id) %>%
		dplyr::filter ( complex_id_a == complex_id_b & complex_id_b == complex_id_c) 
	
	return ( motif_and_protein_complexes)
}

## Function: count_triplet_motifs_in_protein_complexes
## Count the number of triplet motif in which the protein products of all three genes are in the same protein complex
## Input:
## triplet_motifs_list: A table containing the global list of triplet motifs
##	This input table must contain the following columns, 
#   oln_id_a - The Ordered Locus Name of Gene A
#   oln_id_b - The Ordered Locus Name of Gene B
#   oln_id_c - The Ordered Locus Name of Gene C
#   type_ac - The type of interactions between gene A and gene C, one of p, ku, kd, tu, td
#   type_bc - The type of interactions between gene B and gene C, one of p, ku, kd, tu, td  
## protein_complexes: A table containing the list of protein complexes
##	This input table must contain the following columns, 
#   complex_id - A number denoting the protein complex this protein belongs two. A protein can be found in multiple protein complexes
#   oln_id - The Ordered Locus Name corresponding to the protein 
count_triplet_motifs_in_protein_complexes <- function ( triplet_motifs_list, protein_complexes) {
	
	motif_and_protein_complexes <- count_triplet_motifs_in_protein_complexes_helper ( triplet_motifs_list, protein_complexes) 

	# 
	# motif_and_protein_complexes <- left_join ( triplet_motifs_list, protein_complexes, by=c("oln_id_a" = "oln_id") ) %>%
	# 	rename( complex_id_a = complex_id) %>%
	# 	left_join ( protein_complexes, by=c("oln_id_b" = "oln_id") ) %>%
	# 	rename( complex_id_b = complex_id) %>%
	# 	left_join ( protein_complexes, by=c("oln_id_c" = "oln_id") ) %>%
	# 	rename( complex_id_c = complex_id) %>%
	# 	dplyr::filter ( complex_id_a == complex_id_b & complex_id_b == complex_id_c) %>%
	# 
	
		## This following two lines ensure that we don't double count protien complexes. 
		## Sometimes, all three proteins can be found in more than one complexes.
		motif_and_protein_complexes <- motif_and_protein_complexes %>% 
				select(oln_id_a, oln_id_b, oln_id_c, type_ac, type_bc) %>%
				distinct() %>%
				count_triplet_motifs()
	
	return ( motif_and_protein_complexes)
}



#########################################################

count_triplet_motifs_gi_against_edge_list_helper <- function (triplet_motifs_list, edge_list) {
	
	edge_list_updated <- cbind ( edge_list, edge_attribute=rep ( 1, length(edge_list[,1])) )
	
	triplet_motif_and_edge_lists <- left_join ( triplet_motifs_list, edge_list_updated, by=c("oln_id_a" = "oln_id_a", "oln_id_b" = "oln_id_b") ) %>%
		union ( left_join ( triplet_motifs_list, edge_list_updated, by=c("oln_id_a" = "oln_id_b", "oln_id_b" = "oln_id_a") ) ) %>%
		distinct()
	
	a_ge_b <- filter( triplet_motif_and_edge_lists, oln_id_a >= oln_id_b ) %>%
		dplyr::distinct( oln_id_a, oln_id_b, type_ac, type_bc,  edge_attribute) 
	
	b_gt_a <- filter( triplet_motif_and_edge_lists, oln_id_b > oln_id_a ) %>%
		dplyr::distinct( oln_id_a, oln_id_b, type_ac, type_bc,  edge_attribute) %>%
		dplyr::rename ( oln_id_a= oln_id_b, oln_id_b = oln_id_a ) 
	
	edges_in_motifs <- dplyr::union( a_ge_b, b_gt_a)  
	
	return( edges_in_motifs)	
}




## Functipn: count_triplet_motifs_gi_against_edge_list
## Description: Count the number of gene pairs that matches the negatively interaction interaction in each triplet motif 
## Inputs:
## triplet_motifs_list: A table containing the global list of triplet motifs
##	This input table must contain the following columns, 
#   oln_id_a - The Ordered Locus Name of Gene A
#   oln_id_b - The Ordered Locus Name of Gene B
#   oln_id_c - The Ordered Locus Name of Gene C
#   type_ac - The type of interactions between gene A and gene C, one of p, ku, kd, tu, td
#   type_bc - The type of interactions between gene B and gene C, one of p, ku, kd, tu, td  
## edge_list: A table containing the following columns, which represents a relationship between two genes
#   oln_id_a - The Ordered Locus Name of Gene A
#   oln_id_b - The Ordered Locus Name of Gene B

###           This list could also denote any other type of genes e.g. mutants with a particular phenotype etc..)
## Outputs:
## A table containing the following columns:
## type_ac: Type of interaction between gene A and gene C of triplet motif
## type_bc: Type of interaction between gene B and gene C of triplet motif
## count_edge_attributes: Total number of motifs that contains the query set of gene-gene relationships

count_triplet_motifs_gi_against_edge_list <- function (triplet_motifs_list, edge_list) {
	
	edge_list_updated <- cbind ( edge_list, edge_attribute=rep ( 1, length(edge_list[,1])) )
	
	triplet_motif_and_edge_lists <- left_join ( triplet_motifs_list, edge_list_updated, by=c("oln_id_a" = "oln_id_a", "oln_id_b" = "oln_id_b") ) %>%
		union ( left_join ( triplet_motifs_list, edge_list_updated, by=c("oln_id_a" = "oln_id_b", "oln_id_b" = "oln_id_a") ) ) %>%
		distinct()

	a_ge_b <- filter( triplet_motif_and_edge_lists, type_ac >= type_bc ) %>%
		dplyr::distinct( type_ac, type_bc, oln_id_a, oln_id_b, edge_attribute) %>%
		dplyr::group_by( type_ac, type_bc) %>%
		dplyr::summarise( count_edge_attribute=sum(ifelse(is.na(edge_attribute), 0, 1))   )
	
	b_gt_a <- filter( triplet_motif_and_edge_lists, type_bc > type_ac ) %>%
		dplyr::distinct( type_ac, type_bc, oln_id_a, oln_id_b, edge_attribute) %>%
		dplyr::group_by( type_ac, type_bc) %>%
		dplyr::summarise( count_edge_attribute=sum(ifelse(is.na(edge_attribute), 0, 1)) ) %>%
		dplyr::rename ( type_ac= type_bc, type_bc = type_ac ) 
	
	count_edges_in_motifs <- dplyr::union( a_ge_b, b_gt_a)  %>%
		group_by( type_ac, type_bc) %>% 
		summarise(count_edge_attribute=sum(count_edge_attribute) ) %>%
		arrange( type_ac, type_bc)
	
	return( count_edges_in_motifs)	
	
}


## Functipn: count_triplet_motifs_non_gi_edges_against_edge_list
## Description: Count the number of gene pairs that matches the edges separate to the negatively interaction interaction in each triplet motif 
## Inputs:
## triplet_motifs_list: A table containing the global list of triplet motifs
##	This input table must contain the following columns, 
#   oln_id_a - The Ordered Locus Name of Gene A
#   oln_id_b - The Ordered Locus Name of Gene B
#   oln_id_c - The Ordered Locus Name of Gene C
#   type_ac - The type of interactions between gene A and gene C, one of p, ku, kd, tu, td
#   type_bc - The type of interactions between gene B and gene C, one of p, ku, kd, tu, td  
## edge_list: A table containing the following columns, which represents a relationship between two genes
#   oln_id_a - The Ordered Locus Name of Gene A
#   oln_id_b - The Ordered Locus Name of Gene B

###           This list could also denote any other type of genes e.g. mutants with a particular phenotype etc..)
## Outputs:
## A table containing the following columns:
## type_ac: Type of interaction between gene A and gene C of triplet motif
## type_bc: Type of interaction between gene B and gene C of triplet motif
## count_edge_attributes: Total number of motifs that contains the query set of gene-gene relationships
count_triplet_motifs_non_gi_edges_against_edge_list <- function (triplet_motifs_list, edge_list) {
	
	edge_list_updated <- cbind ( edge_list, edge_attribute=rep ( 1, length(edge_list[,1])) )
	
	triplet_motif_and_edge_lists <- left_join ( triplet_motifs_list, edge_list_updated, by=c("oln_id_a" = "oln_id_a", "oln_id_c" = "oln_id_b") ) %>%
		union ( left_join ( triplet_motifs_list, edge_list_updated, by=c("oln_id_a" = "oln_id_b", "oln_id_c" = "oln_id_a") ) ) %>%
		union ( left_join ( triplet_motifs_list, edge_list_updated, by=c("oln_id_b" = "oln_id_a", "oln_id_c" = "oln_id_b") ) ) %>%
		union ( left_join ( triplet_motifs_list, edge_list_updated, by=c("oln_id_b" = "oln_id_b", "oln_id_c" = "oln_id_a") ) ) %>%
		distinct()
	
	## count_a_and_b: If the attribute is found in either gene A or gene B of the triplet motif, I will count this as 1
	## count_c: If the attribute is found in gene C of the triplet motif, I will count this as 1
	
	a_ge_b <- filter( triplet_motif_and_edge_lists, type_ac >= type_bc ) %>%
		dplyr::distinct( type_ac, type_bc, oln_id_a, oln_id_b, edge_attribute) %>%
		dplyr::group_by( type_ac, type_bc) %>%
		dplyr::summarise( count_edge_attribute=sum(ifelse(is.na(edge_attribute),0,1))   )
	
	b_gt_a <- filter( triplet_motif_and_edge_lists, type_bc > type_ac ) %>%
		dplyr::distinct( type_ac, type_bc, oln_id_a, oln_id_b, edge_attribute) %>%
		dplyr::group_by( type_ac, type_bc) %>%
		dplyr::summarise( count_edge_attribute=sum(ifelse(is.na(edge_attribute), 0,1 )) )  %>%
		dplyr::rename ( type_ac= type_bc, type_bc = type_ac) 
	
	count_edges_in_motifs <- dplyr::union( a_ge_b, b_gt_a)  %>%
		group_by( type_ac, type_bc) %>% 
		summarise(count_edge_attribute=sum(count_edge_attribute) ) %>%
		arrange( type_ac, type_bc)
	
	return( count_edges_in_motifs)	
	
}


#########################################################

## Function: count_triplet_motifs_unique_gi_pairs
## Used in identifying repeated use of the negative genetic interactions in the set of all triplet motifs
## Input: 

## triplet_motifs - Table containing the triplet motifs 
##	This input table must contain the following columns, 
#   oln_id_a - The Ordered Locus Name of Gene A
#   oln_id_b - The Ordered Locus Name of Gene B
#   oln_id_c - The Ordered Locus Name of Gene C
#   type_ac - The type of interactions between gene A and gene C, one of p, ku, kd, tu, td
#   type_bc - The type of interactions between gene B and gene C, one of p, ku, kd, tu, td  

## Output:
## A table containing the following columns:
## type_ac: Type of interaction between gene A and gene C of triplet motif
## type_bc: Type of interaction between gene B and gene C of triplet motif
## count: The number of unique negative genetic interactions shared by this type of triplet motif

count_triplet_motifs_unique_gi_pairs <- function (triplet_motifs ) {
	
	a_ge_b <- filter( triplet_motifs, type_ac >= type_bc ) %>%
		group_by ( type_ac, type_bc, oln_id_a, oln_id_b )  %>%
		summarise(counts_level=n()) %>%
		group_by( type_ac, type_bc) %>%
		summarise(counts=n())
	
	b_gt_a <- filter( triplet_motifs, type_bc > type_ac ) %>%
		group_by ( type_ac, type_bc, oln_id_a, oln_id_b )  %>%
		summarise(counts_level=n()) %>%
		group_by( type_ac, type_bc) %>%
		summarise(counts=n()) %>%
		rename ( type_ac= type_bc, type_bc = type_ac) %>%
		select( one_of( c("type_ac", "type_bc", "counts") ))
	
	triplet_motif_counts <- dplyr::union( a_ge_b, b_gt_a)  %>%
		group_by( type_ac, type_bc) %>% 
		summarise(total_count=sum(counts)) %>%
		arrange( type_ac, type_bc)
	
	return( triplet_motif_counts)
}



#########################################################

### Function: get_multiple_edges
# Identify multiple edges in the input network. 

# Inputs: 
# edge_list - A table with the edge list for the network
# source_node_column - The name of the column with the source nodes
# target_node_column - The name of the column with the target nodes

# Output:
# A table with the columns named 'source_node_column' and 'target_node_column'. 
# The table list all of the multiple edges, but with the repeats removed.

get_multiple_edges <- function (edge_list, source_node_column, target_node_column, directed=FALSE, duplicates=FALSE  ) {
  
  edge_list <- edge_list %>% dplyr::select(one_of( c(source_node_column, target_node_column)))
	
  
  multiple_edges <- dplyr::group_by_( edge_list, .dots=c( source_node_column, target_node_column ) ) %>%
  	dplyr::summarise( count=n() ) %>%
  	dplyr::filter(count > 1) %>% 
  	dplyr::select(one_of( c(source_node_column, target_node_column) )) %>% 
  	as.data.frame()
  
  ### Move nodes to one direction only
  if ( directed == FALSE) {
  	
  	rename_condition <- c( target_node_column, source_node_column)
  	names(rename_condition ) <- c(source_node_column, target_node_column)
  	
  	edge_list_reverse <- edge_list  %>% 
  							dplyr::select(one_of( c( target_node_column, source_node_column))) 	%>% 
  							dplyr::rename_( .dots=rename_condition )

  	
  	inner_join_condition <- c( source_node_column, target_node_column)
  	names( inner_join_condition ) <- c( source_node_column, target_node_column)
  	
  	multiple_edges_reverse <-  dplyr::inner_join ( edge_list, edge_list_reverse, by=inner_join_condition) %>%
  							   dplyr::distinct()
  	
  	multiple_edges <- rbind( multiple_edges, multiple_edges_reverse)
  	
  }

  # Get back all the repeated multiple edges
  if ( duplicates == TRUE ) {
  	inner_join_condition <- c( source_node_column, target_node_column)
  	names( inner_join_condition ) <- c( source_node_column, target_node_column)
  	multiple_edges <- dplyr::inner_join ( edge_list, multiple_edges, by=inner_join_condition) %>% 
  					          as.data.frame()
  }
  
  multiple_edges <- multiple_edges %>% arrange_( .dots= c( source_node_column, target_node_column))
  
  return( multiple_edges)
}


### Function: get_multiple_edges
# Identify self-loops in the input network. 

# Inputs: 
# edge_list - A table with the edge list for the network
# source_node_column - The name of the column with the source nodes
# target_node_column - The name of the column with the target nodes

# Output:
# A table with the columns named 'source_node_column' and 'target_node_column'. 
# The table list all of the self-loops in the network.
get_self_loops <- function ( edge_list, source_node_column, target_node_column  ) {
  
  self_loops <-  dplyr::filter_ ( edge_list, .dots= paste ( source_node_column, "==", target_node_column ) ) %>%
    as.data.frame()
  
  return ( self_loops)
}


#########################################################
