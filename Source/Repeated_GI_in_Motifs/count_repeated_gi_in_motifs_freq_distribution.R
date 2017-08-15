### Script: count_triplet_motifs_repeated_gi_edges.R
### Author: Ignatius Pang 
### Date: 27-5-2016
### Description: For each type of triplet motif, count the number of motifs in which multiple triplet motifs of the same type share the same genetic interaction.  
### Repeat the same calculation with randomizated 2000 networks and obtain the bootstrap p-value.
### Network randomization preserves the degree distribution of the network. 

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
		source( "./Source/Common/parameters_file.R")
	}
} else {
	source( "./Source/Common/parameters_file.R")
}

### Local Parameters
if (is_run_locally) {
	results_directory <- "./Results/Bootstrap_p_values_temp/Repeated_GI_in_Motifs_Freq_Dist/"
	create_directory_if_not_exists(results_directory)
}

observed_counts_unique_gi_pairs_in_motifs_file     <- "observed_counts_unique_gi_pairs_in_motifs.tab"

randomized_counts_unique_gi_pairs_in_motifs_file   <- "randomized_counts_unique_gi_pairs_in_motifs.tab"
			  
full_results_counts_unique_gi_pairs_in_motifs_file <- "full_results_counts_unique_gi_pairs_in_motifs.tab"		

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

#### Tallying the number of unique pairs of negative genetic interaction for each type of triplet motif
# observed_counts_unique_gi_pairs_in_motifs <- count_repeated_gi_in_triplet_motifs( triplet_motifs_costanzo)

observed_counts_unique_gi_pairs_in_motifs <- count_repeated_gi_in_triplet_motifs_distribution( triplet_motifs_costanzo)

# observed_counts_unique_gi_pairs_in_motifs <- transpose_count_repeated_gi_in_triplet_motifs_results(observed_counts_unique_gi_pairs_in_motifs, 
# 																 "total_count", "motif_type") 	

write.table( observed_counts_unique_gi_pairs_in_motifs, 
			 file = paste( results_directory, observed_counts_unique_gi_pairs_in_motifs_file, sep=""), 
			 row.names = FALSE)

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
# 																  num_iterations=0.001,
# 																  count_repeated_gi_in_triplet_motifs)

# Need to have mc.set.seed = TRUE and mc.preschedule=FALSE as per user manual on ?mcparallel in the parallel library. See section on Random Numbers.
## "The behaviour with mc.set.seed = TRUE is different only if RNGkind("L'Ecuyer-CMRG") has been selected. Then each time a child is forked it is given the next stream (see nextRNGStream). So if you select that generator, set a seed and call mc.reset.stream just before the first use of mcparallel the results of simulations will be reproducible provided the same tasks are given to the first, second, ... forked process." -- user manual on ?mcparallel in the parallel library, section on Random Numbers

#randomized_counts_orthomcl_paralogs_file <- "randomized_counts_orthomcl_paralogs_in_tri_motifs.tab"
#randomized_counts_sgd_paralogs_file      <- "randomized_counts_sgd_paralogs_in_tri_motifs.tab"

mc.reset.stream() 
randomized_counts_unique_gi_pairs_in_motifs_list <- mclapply ( X=1:number_of_randomized_trials, 
															   FUN= run_one_randomized_trial_compare_with_dataset,  
													  edited_fitered_costanzo_stringent, kinase_network_subset, 
													  sbi_interactome, tf_network, 
													  num_iterations=num_iteration_rewire_network, 
													  count_repeated_gi_in_triplet_motifs_distribution, 
													  mc.set.seed = TRUE , 
													  mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )

names (randomized_counts_unique_gi_pairs_in_motifs_list ) <- 1:length(randomized_counts_unique_gi_pairs_in_motifs_list)

detach(package:dplyr)
library(plyr)
library(dplyr)
randomized_counts_unique_gi_pairs_in_motifs <- ldply( randomized_counts_unique_gi_pairs_in_motifs_list, 
													  as.data.frame ) %>%
												dplyr::rename( group= .id )
detach(package:plyr)
library(dplyr)

write.table ( randomized_counts_unique_gi_pairs_in_motifs, file = paste( results_directory, 
																		 randomized_counts_unique_gi_pairs_in_motifs_file, sep=""), 
			  															 row.names = FALSE)

#########################################################

# # Write the output tables, use the p-value for two-sided test for the OrthoMCL paralogs 
# full_results_counts_unique_gi_pairs_in_motifs <- get_full_results_table(observed_counts_unique_gi_pairs_in_motifs, 
# 																		randomized_counts_unique_gi_pairs_in_motifs, p.value=0.05) 
# 
# write.table ( full_results_counts_unique_gi_pairs_in_motifs, file = paste( results_directory, full_results_counts_unique_gi_pairs_in_motifs_file, sep=""), 
# 			  row.names = TRUE) 

#########################################################

