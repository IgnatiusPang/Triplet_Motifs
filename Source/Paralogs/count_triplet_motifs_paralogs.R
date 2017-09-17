### Script: count_triplet_motifs_paralogs.R
### Author: Ignatius Pang 
### Date: 27-5-2016
### Description: For each type of triplet motif, count the number of motifs in which the edge with the negative genetic interaction 
### consist of a pair of paralogous genes . 
### Repeat the same calculation with randomizated 2000 networks and obtain the bootstrap p-value.
### Randomization conserves the degree distribution of the network. 

#########################################################
# Source location

# psql sbi_triplet_motifs
# cd '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/'
# 

#########################################################
### Parameters 
options <- commandArgs(trailingOnly = TRUE)

### Source the parameters file here 
if ( length(options) != 0  )  { 
	if ( options[1] == 'katana' | options[1] == 'clive')   {
		source( "./parameters_file.R" )

	} else {
		source( "./Common/parameters_file.R")
	}	
} else {
	source( "./Common/parameters_file.R")
}	

### Local Parameters
if (is_run_locally) {

	results_directory <- file.path( results_directory, "Bootstrap_p_values/Paralogs")
	create_directory_if_not_exists(results_directory)
}

observed_counts_orthomcl_paralogs_file	 <- "observed_counts_orthomcl_paralogs_in_tri_motifs.tab"
observed_counts_sgd_paralogs_file	     <- "observed_counts_sgd_paralogs_in_tri_motifs.tab"

randomized_counts_orthomcl_paralogs_file <- "randomized_counts_orthomcl_paralogs_in_tri_motifs.tab"
randomized_counts_sgd_paralogs_file      <- "randomized_counts_sgd_paralogs_in_tri_motifs.tab"
			  
full_results_orthomcl_paralogs_file 	 <- "full_results_table_orthomcl_paralogs_in_tri_motifs.tab"		  
full_results_sgd_paralogs_file  	     <- "full_results_table_sgd_paralogs_in_tri_motifs.tab"

#########################################################
#### Get the actual observed number of triplet motifs

my_join_ac <- c( "query_oln_id_edited" = "oln_id_a")
my_join_bc <- c( "array_oln_id_edited"= "oln_id_a", "oln_id_b" = "oln_id_b" )
my_selected_columns <- c( "query_oln_id_edited", "array_oln_id_edited", 
						  "oln_id_b",  "interaction_type_abbrev.x",
						  "interaction_type_abbrev.y", "genetic_interaction_score", "p_value", "std_dev" )
my_column_names <- c( "oln_id_a", 
					  "oln_id_b",
					  "oln_id_c",
					  "type_ac",
					  "type_bc",
					  "genetic_interaction_score",
					  "p_value",
					  "std_dev" )      

triplet_motifs_costanzo <- form_triplet_motifs(filtered_costanzo_stringent, interactions_combined, interactions_combined,
											   my_join_ac, my_join_bc, my_selected_columns, my_column_names)

#########################################################
#########################################################
#########################################################

## Combine SGD and OrthoMCL datasets
# sgd_paralogs <- filter( sgd_paralogs, oln_id_a < oln_id_b )  %>%
# 	dplyr::select ( one_of ( c("oln_id_b", "oln_id_a")) ) %>%
# 	dplyr::rename(  oln_id_c=oln_id_a ) %>% 
# 	dplyr::rename( oln_id_a=oln_id_b) %>%
# 	dplyr::rename ( oln_id_b=oln_id_c) %>%
# 	dplyr::union (  filter( sgd_paralogs, oln_id_a > oln_id_b ) )  %>%
# 	dplyr::distinct()
# 
# orthomcl_paralogs <- filter( orthomcl_paralogs, oln_id_a < oln_id_b )  %>%
# 	dplyr::select ( one_of ( c("oln_id_b", "oln_id_a")) ) %>%
# 	dplyr::rename(  oln_id_c=oln_id_a ) %>% 
# 	dplyr::rename( oln_id_a=oln_id_b) %>%
# 	dplyr::rename ( oln_id_b=oln_id_c) %>%
# 	dplyr::union (  filter( orthomcl_paralogs, oln_id_a > oln_id_b ) ) %>% 
# 	dplyr::distinct()
# 

orthomcl_paralogs <- clean_graph_table  (orthomcl_paralogs,  "oln_id_a", "oln_id_b", directed=FALSE ) 

sgd_paralogs 	  <- clean_graph_table  (sgd_paralogs,  "oln_id_a", "oln_id_b", directed=FALSE ) 

orthomcl_paralogs <- dplyr::distinct ( dplyr::union ( orthomcl_paralogs, sgd_paralogs) )

orthomcl_paralogs <- clean_graph_table  (orthomcl_paralogs,  "oln_id_a", "oln_id_b", directed=FALSE ) 

##### Analysis to help count the number of proteins involved in paralogs 
orthomcl_paralogs_table <- as.data.frame ( orthomcl_paralogs )
length( unique ( c( orthomcl_paralogs_table[, "oln_id_a"], orthomcl_paralogs_table[, "oln_id_b"] ) ) )
# 1542  proteins in total 

#### How many of the pairs of paralogs / ohnologs shared negative genetic interactions?


gi_and_paralog_table_a <- inner_join ( orthomcl_paralogs_table, filtered_costanzo_stringent_2016, 
									   								by= c("oln_id_a" = "query_oln_id_edited", 
																		  "oln_id_b" = "array_oln_id_edited" )) %>%
	dplyr::select(one_of (c('oln_id_a', 'oln_id_b')) )

gi_and_paralog_table_b <- inner_join ( orthomcl_paralogs_table, filtered_costanzo_stringent_2016 , 
									   								by= c("oln_id_b" = "query_oln_id_edited", 
																		  "oln_id_a" = "array_oln_id_edited" ))  %>%
	dplyr::select(one_of (c('oln_id_a', 'oln_id_b')) ) 

result_table <- dplyr::union ( gi_and_paralog_table_a, gi_and_paralog_table_b)

##### 

# length ( sgd_paralogs[,1]) 
# length ( orthomcl_paralogs[,1]) 
# 
# count ( sgd_paralogs[,1])
# count ( orthomcl_paralogs[,1])

#### Analysis on OrthoMCL list of paralogs 

# Test the edge that is involved in negative genetic interactions
count_motifs_with_orthomcl_paralogs <- count_triplet_motifs_gi_against_edge_list( triplet_motifs_costanzo, orthomcl_paralogs)

count_motifs_with_orthomcl_paralogs <- transpose_count_triplet_motifs_results( count_motifs_with_orthomcl_paralogs, 
																"type_ac", "type_bc", "count_edge_attribute", "motif_type") 	

write.table( count_motifs_with_orthomcl_paralogs, file = file.path( results_directory, observed_counts_orthomcl_paralogs_file), 
			 row.names = FALSE)


### Helper function to get statistics for paper.
## How many triplet motifs are associated with pairs of genes that are paralogs?
list_of_motifs_with_paralogs <- count_triplet_motifs_gi_against_edge_list_helper (triplet_motifs_costanzo, orthomcl_paralogs)
list_of_motifs_with_paralogs


# How many unique pairs of paralogs are associated with triplet motifs?
count_of_motifs_with_paralogs <- list_of_motifs_with_paralogs %>%
								 dplyr::filter(edge_attribute == 1) %>%
								 dplyr::select( one_of( c("oln_id_a", "oln_id_b") ) ) %>%
								 dplyr::distinct( oln_id_a, oln_id_b ) %>%
								 as.data.frame() %>%
								 dplyr::count()
count_of_motifs_with_paralogs


#########################################################
#########################################################
#########################################################
#########################################################

### Randomized triplet motifs analysis 10 times, example
### This uses the mclappy function to spead the calculations onto different cores. It will only work on Linux machines and not Windows.

edited_fitered_costanzo_stringent <- filtered_costanzo_stringent

colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "query_oln_id_edited"] <- "oln_id_a"
colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "array_oln_id_edited"] <- "oln_id_b"

### Run one trial to see whether it works or not 
# result_one_trial <- run_one_randomized_trial_compare_with_dataset(1, edited_fitered_costanzo_stringent, 
# 																  kinase_network_subset, sbi_interactome, tf_network, 
# 																  num_iterations=100,
# 																  count_triplet_motifs_gi_against_edge_list, 
#													  			  min_num_neighbours=10, max_num_tries=10, mode="dd",
# 																  orthomcl_paralogs )

# Need to have mc.set.seed = TRUE and mc.preschedule=FALSE as per user manual on ?mcparallel in the parallel library. See section on Random Numbers.
## "The behaviour with mc.set.seed = TRUE is different only if RNGkind("L'Ecuyer-CMRG") has been selected. Then each time a child is forked it is given the next stream (see nextRNGStream). So if you select that generator, set a seed and call mc.reset.stream just before the first use of mcparallel the results of simulations will be reproducible provided the same tasks are given to the first, second, ... forked process." -- user manual on ?mcparallel in the parallel library, section on Random Numbers

#randomized_counts_orthomcl_paralogs_file <- "randomized_counts_orthomcl_paralogs_in_tri_motifs.tab"
#randomized_counts_sgd_paralogs_file      <- "randomized_counts_sgd_paralogs_in_tri_motifs.tab"

mc.reset.stream() 
randomized_counts_orthomcl_paralogs_list <- mclapply ( X=1:number_of_randomized_trials, FUN= run_one_randomized_trial_compare_with_dataset,   
													  edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
													  num_iterations=num_iteration_rewire_network, 
													  count_triplet_motifs_gi_against_edge_list, 
													  min_num_neighbours=10, max_num_tries=10, mode="dd",													  
													  orthomcl_paralogs,
													  mc.set.seed = TRUE , 
													  mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )
											  
randomized_counts_orthomcl_paralogs_table <- concat_motif_counts_list_into_table( randomized_counts_orthomcl_paralogs_list) 

write.table ( randomized_counts_orthomcl_paralogs_table, file = paste( results_directory, randomized_counts_orthomcl_paralogs_file), 
			  row.names = FALSE)

# ###
# randomized_counts_sgd_paralogs_list <- mclapply ( X=1:number_of_randomized_trials, FUN= run_one_randomized_trial_compare_with_dataset,   
														  # edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
												  		  # num_iterations=num_iteration_rewire_network, 
														  # count_triplet_motifs_gi_against_edge_list, 
														  # min_num_neighbours=10, max_num_tries=10, mode="dd",													  														  
														  # sgd_paralogs,
														  # mc.set.seed = TRUE, 
														  # mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )

# randomized_counts_sgd_paralogs_table <- concat_motif_counts_list_into_table( randomized_counts_sgd_paralogs_list) 

# write.table ( randomized_counts_sgd_paralogs_table, file = paste( results_directory, randomized_counts_sgd_paralogs_file), 
# row.names = FALSE)

#########################################################

# Write the output tables, use the p-value for two-sided test for the OrthoMCL paralogs 
full_results_orthomcl_paralogs <- get_full_results_table(count_motifs_with_orthomcl_paralogs, 
														 randomized_counts_orthomcl_paralogs_table, p.value=0.05) 

write.table ( full_results_orthomcl_paralogs, file = paste( results_directory, full_results_orthomcl_paralogs_file), 
			  row.names = TRUE) 

# # Write the output tables, use the p-value for two-sided test for the SGD paralogs 
# full_results_sgd_paralogs <- get_full_results_table(count_motifs_with_sgd_paralogs, 
													# randomized_counts_sgd_paralogs_table, p.value=0.05) 

# write.table ( full_results_sgd_paralogs, file = paste( results_directory, full_results_sgd_paralogs_file), 
			  # row.names = TRUE) 

#########################################################
