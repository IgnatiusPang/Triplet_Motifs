### Script: count_triplet_motifs_negative_interactions.R
### Author: Ignatius Pang 
### Date: 20-5-2016
### Description: Count the observed number of triplet motifs. Randomize the respective networks and count the number triplet motifs, repeat randomization 1000 times.
### Obtain the bootstrap p-value

##########################################################
# Source location

# psql sbi_triplet_motifs
# cd '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/'
# 

##########################################################
### Global Parameters 

options <- commandArgs(trailingOnly = TRUE)

### Source the header file here 

if ( length(options) != 0 
	 & ( options[1] == 'katana' | options[1] == 'clive')  ) {

			source( "./parameters_file.R" )

} else {
	source( "./Common/parameters_file.R")
}

### Local Parameters

output_observed_triplet_motif_counts  <- "results_count_triplet_motifs_costanzo.tab" # "results_count_triplet_motifs_costanzo.tab"
randomized_results_table_file         <- "randomized_results_table.tab"
output_enrichment_p_value             <- "enrichment_p_value_fixed.tab"
output_depletion_p_value              <- "depletion_p_value_fixed.tab"
output_full_results_table             <- "full_results_table.tab"

if (is_run_locally) {

	results_directory <- file.path( results_directory, "Bootstrap_p_values_temp/Negative_Genetic_Interactions")
	
	create_directory_if_not_exists(results_directory)
}


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

results_count_triplet_motifs_costanzo <- count_triplet_motifs_custom(triplet_motifs_costanzo, "type_ac", "type_bc", "my_total_count")

results_count_triplet_motifs_costanzo <- transpose_count_triplet_motifs_results  (results_count_triplet_motifs_costanzo,  "type_ac", "type_bc", "my_total_count", "motif_type") 
	
write.table(results_count_triplet_motifs_costanzo , file =  file.path( results_directory, output_observed_triplet_motif_counts), 
			row.names = FALSE)

if (is_run_locally) {
	
	if ( use_costanzo_2010_dataset == TRUE) {
		# Backup all the triplet motifs in the data directory
		saveRDS( triplet_motifs_costanzo, file=file.path(list_of_triplets_directory, "triplet_motifs_costanzo_2010.Rdata") ) 
		write.table( triplet_motifs_costanzo, file=file.path(list_of_triplets_directory, "triplet_motifs_costanzo_2010.txt"), row.names=FALSE ) 
		
	} else {
		saveRDS( triplet_motifs_costanzo, file=file.path(list_of_triplets_directory, "triplet_motifs_costanzo_2016.Rdata") ) 
		write.table( triplet_motifs_costanzo, file=file.path(list_of_triplets_directory, "triplet_motifs_costanzo_2016.txt"), row.names=FALSE ) 
		
	}
}

#########################################################

### Randomized triplet motifs analysis 10 times, example
### This uses the mclappy function to spead the calculations onto different cores. It will only work on Linux machines and not Windows.

edited_fitered_costanzo_stringent <- filtered_costanzo_stringent

colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "query_oln_id_edited"] <- "oln_id_a"
colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "array_oln_id_edited"] <- "oln_id_b"


### Run one trial to see whether it works or not 
# result_one_trial <- run_one_randomized_trial_compare_with_dataset(1, edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network,
#  									num_iterations=num_iteration_rewire_network, count_triplet_motifs)


# Need to have mc.set.seed = TRUE and mc.preschedule=FALSE as per user manual on ?mcparallel in the parallel library. See section on Random Numbers.
## "The behaviour with mc.set.seed = TRUE is different only if RNGkind("L'Ecuyer-CMRG") has been selected. Then each time a child is forked it is given the next stream (see nextRNGStream). 
## So if you select that generator, set a seed and call mc.reset.stream just before the first use of mcparallel the results of simulations will be reproducible provided the same tasks are given to the first, second, ... forked process." 
##  -- user manual on ?mcparallel in the parallel library, section on Random Numbers

mc.reset.stream() 
list_of_randomized_triplet_motif_counts <- mclapply ( X=1:number_of_randomized_trials, FUN= run_one_randomized_trial_compare_with_dataset,   
													  edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
													  num_iterations=num_iteration_rewire_network, 
													  count_triplet_motifs,
													  mc.set.seed = TRUE , 
													  mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )

 randomized_results_table <- concat_motif_counts_list_into_table( list_of_randomized_triplet_motif_counts) 

write.table ( randomized_results_table, file = file.path( results_directory, randomized_results_table_file), 
			  row.names = FALSE)

#########################################################

## Calculate the p-values from the enrichment and depletion analyses

bootstrap_p_values <- calculate_boostrap_p_values( results_count_triplet_motifs_costanzo, randomized_results_table ) 

write.table ( bootstrap_p_values$enriched, file = file.path( results_directory, output_enrichment_p_value), 
			  row.names = TRUE) 

write.table ( bootstrap_p_values$depleted, file = file.path( results_directory, output_depletion_p_value), 
			  row.names = TRUE) 

#########################################################

# Write the output tables, use the p-value for two-sided test
final_results_table <- get_full_results_table(results_count_triplet_motifs_costanzo, randomized_results_table, p.value=0.05) 

write.table ( final_results_table, file = file.path( results_directory, output_full_results_table), 
			  row.names = TRUE) 

#########################################################

# results_directory <- "/media/z3371724/PostDoc/2016/Triplet_Motifs/Results/Bootstrap_p_values/Negative_Genetic_Interactions/"
# randomized_results_table <- read.table (  file = file.path( results_directory, randomized_results_table_file), header = TRUE)
# results_count_triplet_motifs_costanzo <- read.table(file =  file.path( results_directory, output_observed_triplet_motif_counts),  header = TRUE)
# 
# 
  
  
# SNF1 / YDR477W Overview, involved in tuku interactions. i.e. it is a kinase but also have transcription factor regulatory functions.
