### Script: compare_triplet_motif_gi_with_other_networks.R
### Author: Ignatius Pang 
### Date: 22-8-2016
### Description: Compare the genetic interaction edge in triplet motifs with other types of interactions. 

#########################################################
### Parameters 
### Global Parameters 

options <- commandArgs(trailingOnly = TRUE)

debug <- TRUE

### Source the header file here 

if ( length(options) != 0 
	 & ( options[1] == 'katana' | options[1] == 'clive')  ) {
	source( "./parameters_file.R" )

} else {
	source( "./Common/parameters_file.R")
}

### Local Parameters
observed_results_table_file   <- "observed_results_table_file.tab"
randomized_results_table_file <- "randomized_results_table.tab"
output_full_results_table     <- "full_results_table.tab"

if (is_run_locally) {
	
	results_directory <- file.path( results_directory, "Bootstrap_p_values_temp/Analyze_GI_edges" ) 
	
	create_directory_if_not_exists(final_results_directory)
	
}

#########################################################

#### Get the actual observed number of triplet motifs


#########################################################

fix_edge_attribute_direction <- function (interaction_type_abbrev) {
	
	if ( is.na( interaction_type_abbrev ) ) {
		return (NA)
	} 	
	
	if (interaction_type_abbrev== "tu" ) {
		return ('td')
	} else if ( interaction_type_abbrev== "ku") {
		return ('kd')
	} else if ( interaction_type_abbrev== "td") {
		return ('tu')
	} else if ( interaction_type_abbrev== "kd") {
		return ('ku')
	} else if ( interaction_type_abbrev== "p") {
		return ('p')
	} 
	else {
		stop ('interaction type not supported')
	}
}


#########################################################

# negative_genetic_interaction <- filtered_costanzo_stringent
# oln_id_a <- 'oln_id_a'
# oln_id_b <- 'oln_id_b'
# type_ac <- 'type_ac'
# type_bc <- 'type_bc'
# interaction_type_abbrev <- 'interaction_type_abbrev'
# my_counts <- 'count'
# total_count <- 'my_total_count'

count_triplet_motifs_compare_with_other_networks_second_step <- function (negative_genetic_interaction, interactions_combined, oln_id_a, oln_id_b, type_ac, type_bc, interaction_type_abbrev, my_counts, total_count ) {

	colnames( negative_genetic_interaction)[colnames( negative_genetic_interaction) == "oln_id_a"] <-  "query_oln_id_edited"
	colnames( negative_genetic_interaction)[colnames( negative_genetic_interaction) == "oln_id_b"] <-  "array_oln_id_edited"

	my_join_ac <- c( "query_oln_id_edited" = "oln_id_a")
	my_join_bc <- c( "array_oln_id_edited"= "oln_id_a", "oln_id_b" = "oln_id_b" )
	my_selected_columns <- c("query_oln_id_edited", "array_oln_id_edited", 
							 "oln_id_b",  "interaction_type_abbrev.x",
							 "interaction_type_abbrev.y", "genetic_interaction_score")
	my_column_names <- c( "oln_id_a", 
						  "oln_id_b",
						  "oln_id_c",
						  "type_ac",
						  "type_bc",
						  "genetic_interaction_score" )      
	
	triplet_motifs <- form_triplet_motifs(negative_genetic_interaction, interactions_combined, interactions_combined,
												   my_join_ac, my_join_bc, my_selected_columns, my_column_names)

	table_join_by_names <- c( oln_id_a, oln_id_b)
	names( table_join_by_names ) <- c( oln_id_a, oln_id_b)
	
	compare_gi_with_other_networks <- dplyr::left_join(triplet_motifs, interactions_combined,
			  								     		by=table_join_by_names) 
	
	
	if (debug == TRUE) {
		
		triplet_motifs %>%
			select (  one_of ( c('oln_id_a', 'oln_id_b'))) %>%
			filter ( oln_id_a > oln_id_b) %>%
			distinct() %>%
			count() %>% 
			print ( )
		
		filter ( compare_gi_with_other_networks, is.na(interaction_type_abbrev)== TRUE ) %>%
			select (  one_of ( c('oln_id_a', 'oln_id_b'))) %>%
			filter ( oln_id_a > oln_id_b) %>%
			distinct() %>%
			count() %>% 
			print ( )
		
		filter ( compare_gi_with_other_networks, is.na(interaction_type_abbrev)== FALSE ) %>%
		select (  one_of ( c('oln_id_a', 'oln_id_b'))) %>%
			filter ( oln_id_a > oln_id_b) %>%
			distinct() %>%
			count() %>% 
			print ( )
	}

	
	compare_gi_with_other_networks <- compare_gi_with_other_networks %>% 
										dplyr::group_by_ ( .dots=lapply ( c(type_ac, type_bc, interaction_type_abbrev), as.name) ) %>%
										dplyr::summarise(count=n())
		
	compare_gi_with_other_networks <- count_triplet_motifs_compare_with_other_networks( compare_gi_with_other_networks, 
													  type_ac, type_bc, interaction_type_abbrev, my_counts, total_count)
	
	compare_gi_with_other_networks <- as.data.frame( compare_gi_with_other_networks)
	
	compare_gi_with_other_networks[is.na(compare_gi_with_other_networks[, interaction_type_abbrev]), interaction_type_abbrev] <- "NA"
	
	compare_gi_with_other_networks_final_results <- compare_gi_with_other_networks %>% 
									  				tidyr::spread( interaction_type_abbrev, my_total_count ) %>%
													as.data.frame()
	
	return(compare_gi_with_other_networks_final_results )
}

#########################################################
# triplet_motifs <- compare_gi_with_other_networks 


count_triplet_motifs_compare_with_other_networks <- function (triplet_motifs, type_ac, type_bc, interaction_type_abbrev, my_counts, total_count ) {
	
	a_ge_b <- dplyr::filter_( triplet_motifs, paste(type_ac, ">=", type_bc) ) %>%
			  dplyr::group_by_( .dots=lapply ( c(type_ac, type_bc, interaction_type_abbrev), as.name)  ) %>%
			  dplyr::summarise_( .dots= list(counts=interp( ~sum(x), x=as.name(my_counts))) )
	
	columns_to_select <- c(type_ac, type_bc, interaction_type_abbrev, "counts")
	
	b_gt_a <- dplyr::filter_( triplet_motifs, paste(type_bc, ">", type_ac) ) %>%
			  dplyr::group_by_( .dots=lapply ( c(type_ac, type_bc, interaction_type_abbrev), as.name) ) %>%
			  dplyr::summarise_( .dots= list(counts=interp( ~sum(x), x=as.name(my_counts))) ) %>%
			  dplyr::rename_( .dots=setNames(list(type_ac, type_bc, interaction_type_abbrev), 
			  							   c(type_bc, type_ac, interaction_type_abbrev) ) ) %>%
		      dplyr::select( one_of( columns_to_select ))
	
	b_gt_a <- as.data.frame( b_gt_a )
	
	b_gt_a[,interaction_type_abbrev] <- sapply(b_gt_a[,interaction_type_abbrev],
											   fix_edge_attribute_direction)
	
	b_gt_a <- tibble::as_tibble(b_gt_a)
	
	triplet_motif_counts <- dplyr::bind_rows( a_ge_b, b_gt_a)  %>%
		dplyr::group_by_( .dots=lapply ( c(type_ac, type_bc, interaction_type_abbrev), as.name) ) %>%	
		dplyr::summarise( temp_counts=sum(counts))  %>%
		dplyr::rename_( .dots=setNames(list("temp_counts"), c(total_count) ) ) %>%
		dplyr::arrange_( .dots=c(type_ac, type_bc) )
	
	# as.data.frame( a_ge_b ) %>% filter ( type_ac %in% c('td', 'p') & type_bc %in% c('td', 'p') & type_ac != type_bc)
	# as.data.frame( b_gt_a ) %>% filter ( type_ac %in% c('td', 'p') & type_bc %in% c('td', 'p') & type_ac != type_bc)
	# as.data.frame( triplet_motifs ) %>% filter ( type_ac %in% c('td', 'p') & type_bc %in% c('td', 'p') & type_ac != type_bc)
	# as.data.frame( triplet_motif_counts ) %>% filter ( type_ac %in% c('td', 'p') & type_bc %in% c('td', 'p') & type_ac != type_bc)
	
	return( triplet_motif_counts)
}

#########################################################

compare_gi_with_other_networks_final_results <- count_triplet_motifs_compare_with_other_networks_second_step( filtered_costanzo_stringent, 
																											  interactions_combined, 
																											  "oln_id_a", "oln_id_b",
																											  "type_ac", "type_bc", 
																											  "interaction_type_abbrev", 
																											  "count", "my_total_count")

write.table ( compare_gi_with_other_networks_final_results, file=file.path(results_directory, 
																	   observed_results_table_file), row.names = FALSE)

#########################################################

### Randomized triplet motifs analysis 10 times, example
### This uses the mclappy function to spead the calculations onto different cores. It will only work on Linux machines and not Windows.

edited_fitered_costanzo_stringent <- filtered_costanzo_stringent

colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "query_oln_id_edited"] <- "oln_id_a"
colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "array_oln_id_edited"] <- "oln_id_b"


### Run one trial to see whether it works or not 
# result_one_trial <- run_one_randomized_trial_compare_with_dataset(1, edited_fitered_costanzo_stringent, kinase_network_subset, 
# 																  sbi_interactome, tf_network,
# 																  num_iterations=num_iteration_rewire_network, 
# 																  count_triplet_motifs_compare_with_other_networks_second_step, 
# 																  min_num_neighbours=10, max_num_tries=10, mode="dd", 
# 																  use_rewired_interactions=TRUE, "oln_id_a", "oln_id_b",
# 																  "type_ac", "type_bc", "interaction_type_abbrev", 
# 																  "count", "my_total_count" )

#########################################################

mc.reset.stream() 
list_of_randomized_triplet_motif_counts <- mclapply ( X=1:number_of_randomized_trials, FUN= run_one_randomized_trial_compare_with_dataset,   
													  edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
													  num_iterations=num_iteration_rewire_network, 
													  
													  count_triplet_motifs_compare_with_other_networks_second_step, 
													  min_num_neighbours=10, max_num_tries=10, mode="dd", 
													  use_rewired_interactions=TRUE, "oln_id_a", "oln_id_b",
													  "type_ac", "type_bc", "interaction_type_abbrev", 
													  "count", "my_total_count", 
													  
													  mc.set.seed = TRUE , 
													  mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )



#### Concatenate the results together 
# rbind_list_elements_recursively <- function ( x) {
# 	
# 	if (length(x)== 1 ) {
# 		return ( x[[1]])
# 	} else {
# 	
# 		return ( rbind ( x[[1]],  rbind_list_elements_recursively(x[2:length(x)]) ))
# 	}
# }



### Concatenate the results together 
rbind_list_elements_recursively <- function ( x) {
=======
compare_gi_with_other_networks <- compare_gi_with_other_networks %>% 
									dplyr::group_by ( type_ac, type_bc, interaction_type_abbrev ) %>%
									dplyr::summarise(count=n())
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
	
	if (length(x)== 1 ) {
		return ( x[[1]])
	} else {
		
		
		recursive_result <- rbind_list_elements_recursively(x[2:length(x)])
		
		curr_row_colnames <- colnames( x[[1]] )
		prev_row_colnames <- colnames( recursive_result )
		
		if ( is.null( curr_row_colnames ) ) {
			curr_row_colnames <- c()	
		}
		
		if ( is.null( prev_row_colnames ) ) {
			prev_row_colnames <- c()	
		}
		
		in_prev_not_curr <- setdiff( prev_row_colnames, curr_row_colnames )
		in_curr_not_prev <- setdiff( curr_row_colnames, prev_row_colnames )
		
		for ( i in in_prev_not_curr) {
			
			x[[1]] <- cbind( x[[1]],  rep(0, 1) )
			
		}
		
		colnames( x[[1]]) <- c( curr_row_colnames, in_prev_not_curr)
		
		
		recursive_result_row_length <- dim  ( recursive_result )[1]
		
		if ( is.null(recursive_result_row_length) ) {
			recursive_result_row_length <- 1
			
		}
		
		for ( i in in_curr_not_prev) {
			
			recursive_result <- cbind( recursive_result,  rep(0, length(recursive_result_row_length )) )
			
		}
		
		colnames( recursive_result) <- c( prev_row_colnames, in_curr_not_prev)
		
		
		recursive_result <- recursive_result [, colnames( x[[1]]) ]
		
		
		return ( rbind (x[[1]] ,   recursive_result))
	}
}



randomized_results_table <- rbind_list_elements_recursively(list_of_randomized_triplet_motif_counts)

write.table ( randomized_results_table, file = file.path( results_directory, randomized_results_table_file), 
			  row.names = FALSE)

#########################################################

