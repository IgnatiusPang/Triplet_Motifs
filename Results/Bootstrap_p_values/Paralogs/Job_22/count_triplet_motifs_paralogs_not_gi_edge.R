### Script: count_triplet_motifs_paralogs_not_gi.R
### Author: Ignatius Pang 
### Date: 27-5-2016
### Description: For each type of triplet motif, count the number of motifs in which the edge with the negative genetic interaction 
### consist of a pair of paralogous genes . 
### Repeat the same calculation with randomizated 2000 networks and obtain the bootstrap p-value.
### Randomization conserves the degree distribution of the network. 

library(igraph)
library(dplyr)
library(lazyeval)
library(parallel)

sessionInfo()

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
		source( "/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Common/parameters_file.R")
	}
}

### Local Parameters
if (is_run_locally) {
	results_directory <- paste( base_directory, "Results/Bootstrap_p_values_temp/Paralogs/",  sep="")
}

observed_counts_orthomcl_paralogs_not_gi_file 	<- "observed_counts_orthomcl_paralogs_not_gi_in_tri_motifs.tab"
observed_counts_sgd_paralogs_not_gi_file	    <- "observed_counts_sgd_paralogs_not_gi_in_tri_motifs.tab"

randomized_counts_orthomcl_paralogs_not_gi_file <- "randomized_counts_orthomcl_paralogs_not_gi_in_tri_motifs.tab"
randomized_counts_sgd_paralogs_not_gi_file      <- "randomized_counts_sgd_paralogs_not_gi_in_tri_motifs.tab"
			  
full_results_orthomcl_paralogs_not_gi_file    	<- "full_results_table_orthomcl_paralogs_not_gi_in_tri_motifs.tab"		  
full_results_sgd_paralogs_not_gi_file  		    <- "full_results_table_sgd_paralogs_not_gi_in_tri_motifs.tab"

#########################################################

#### Get the actual observed number of triplet motifs

my_join_ac <- c( "query_oln_id_edited" = "oln_id_a")
my_join_bc <- c( "array_oln_id_edited"= "oln_id_a", "oln_id_b" = "oln_id_b" )
my_selected_columns <- c("query_oln_id_edited", "array_oln_id_edited", 
						 "oln_id_b",  "interaction_type_abbrev.x",
						 "interaction_type_abbrev.y", "genetic_interaction_score", "p_value", "std_dev")
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

## Combine SGD ohnologs and OrthoMCL datasets
orthomcl_paralogs <- dplyr::distinct ( dplyr::union ( orthomcl_paralogs, sgd_paralogs) )

orthomcl_paralogs <- clean_graph_table  (orthomcl_paralogs,  "oln_id_a", "oln_id_b", directed=FALSE ) 

#### Analysis on OrthoMCL list of paralogs 


# Test the edges that is not the negative genetic interaction

count_motifs_with_orthomcl_paralogs_not_gi <- count_triplet_motifs_non_gi_edges_against_edge_list (triplet_motifs_costanzo, orthomcl_paralogs) 
	
count_motifs_with_orthomcl_paralogs_not_gi <- transpose_count_triplet_motifs_results(count_motifs_with_orthomcl_paralogs_not_gi, 
																				  "type_ac", "type_bc", "count_edge_attribute", "motif_type") 	


write.table(count_motifs_with_orthomcl_paralogs_not_gi , file =  paste( results_directory, observed_counts_orthomcl_paralogs_not_gi_file, sep=""), 
			row.names = FALSE)


#########################################################
#########################################################
#########################################################
#########################################################

### Randomized triplet motifs analysis 10 times, example
### This uses the mclappy function to spead the calculations onto different cores. It will only work on Linux machines and not Windows.

edited_filtered_costanzo_stringent <- filtered_costanzo_stringent

colnames( edited_filtered_costanzo_stringent)[colnames( edited_filtered_costanzo_stringent) == "query_oln_id_edited"] <- "oln_id_a"
colnames( edited_filtered_costanzo_stringent)[colnames( edited_filtered_costanzo_stringent) == "array_oln_id_edited"] <- "oln_id_b"

### Run one trial to see whether it works or not 
 result_one_trial <- run_one_randomized_trial_compare_with_dataset(1, edited_filtered_costanzo_stringent, 
 																  kinase_network_subset, sbi_interactome, tf_network, 
 																  num_iterations=100,
 																  count_triplet_motifs_non_gi_edges_against_edge_list, 
													   			  min_num_neighbours=10, max_num_tries=10, mode="dd",
 																  orthomcl_paralogs )

# Need to have mc.set.seed = TRUE and mc.preschedule=FALSE as per user manual on ?mcparallel in the parallel library. See section on Random Numbers.
## "The behaviour with mc.set.seed = TRUE is different only if RNGkind("L'Ecuyer-CMRG") has been selected. Then each time a child is forked it is given the next stream (see nextRNGStream). So if you select that generator, set a seed and call mc.reset.stream just before the first use of mcparallel the results of simulations will be reproducible provided the same tasks are given to the first, second, ... forked process." -- user manual on ?mcparallel in the parallel library, section on Random Numbers



mc.reset.stream() 
randomized_counts_orthomcl_paralogs_not_gi_list <- mclapply ( X=1:number_of_randomized_trials, FUN= run_one_randomized_trial_compare_with_dataset,   
													  edited_filtered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
													  num_iterations=num_iteration_rewire_network, 
													  count_triplet_motifs_non_gi_edges_against_edge_list, 
													  min_num_neighbours=10, max_num_tries=10, mode="dd",
													  orthomcl_paralogs,
													  mc.set.seed = TRUE , 
													  mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )

													  

randomized_counts_orthomcl_paralogs_not_gi_table <- concat_motif_counts_list_into_table( randomized_counts_orthomcl_paralogs_not_gi_list) 

write.table ( randomized_counts_orthomcl_paralogs_not_gi_table, file = paste( results_directory, randomized_counts_orthomcl_paralogs_not_gi_file, sep=""), 
			  row.names = FALSE)




#########################################################

# Write the output tables, use the p-value for two-sided test for the OrthoMCL paralogs 
full_results_orthomcl_paralogs_not_gi <- get_full_results_table(count_motifs_with_orthomcl_paralogs_not_gi, 
														 randomized_counts_orthomcl_paralogs_not_gi_table, p.value=0.05) 

write.table ( full_results_orthomcl_paralogs_not_gi, file = paste( results_directory, full_results_orthomcl_paralogs_not_gi_file, sep=""), 
			  row.names = TRUE) 



#########################################################
