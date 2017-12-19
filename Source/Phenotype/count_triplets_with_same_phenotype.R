### Script: count_triplet_motifs_negative_gi_essential_genes.R
### Author: Ignatius Pang 
### Date: 27-5-2016
### Description: Count the observed number of triplet motifs in which the protein products of all three genes are found in the same protein complex. 
### Repeat the same calculation with randomizated 2000 networks and obtain the bootstrap p-value.
### Randomization conserves the degree distribution of the network. 

#########################################################
# Source location

# psql sbi_triplet_motifs
# cd '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/'
# 

#########################################################
### Global Parameters 

options <- commandArgs(trailingOnly = TRUE)

### Source the header file here 

if ( length(options) != 0 
	 & ( options[1] == 'katana' | options[1] == 'clive')  ) {
	
	source( "./parameters_file.R" )
	
} else {
	source( "./Common/parameters_file.R")
}

supplementary_data_directory <- file.path(results_directory, "Supplementary_Files" ) 
create_directory_if_not_exists(supplementary_data_directory)

### Local Parameters
if (is_run_locally) {

	results_directory <- file.path ( results_directory, "Bootstrap_p_values_temp/Phenotype") 
	create_directory_if_not_exists( results_directory)
}

observed_counts_phenotypes_file    <- paste( "observed_counts_phenotypes_in_tri_motifs.tab", sep="")
randomized_counts_phenotypes_file  <- paste( "randomized_counts_phenotypes_in_tri_motifs.tab", sep="")
full_results_phenotypes_file  	 <- paste( "full_results_table_phenotypes_in_tri_motifs.tab", sep="")

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

##################################################################################################################################################

if (   length(options) == 0 )  {
	
motif_and_phenotypes_detail <- count_triplet_motifs_phenotypes_helper (  triplet_motifs_costanzo, filtered_phenotype_data_detailed, 
																		 is_keep_phenotype_list = TRUE) 

motif_and_phenotypes_detail <- motif_and_phenotypes_detail %>% dplyr::mutate( phenotype_intersect_abc = 
																			  	unlist(map( phenotype_intersect_abc, function(x) { 
																			  		paste( as.character(x), collapse="; ")
																			  		}))  )

motif_and_phenotypes_detail <- motif_and_phenotypes_detail %>% dplyr::mutate( num_phenotype_shared = 
																			  	unlist(map( num_phenotype_shared, function(x) { 
																			  		paste( as.character(x), collapse="; ")
																			  	}))  )

write.table (motif_and_phenotypes_detail, 
			 file=file.path( supplementary_data_directory, "motif_and_phenotypes_detail.txt"), sep="\t", quote = FALSE, row.names = FALSE )

}

##################################################################################################################################################

motif_and_phenotypes <- count_triplet_motifs_phenotypes (  triplet_motifs_costanzo, filtered_phenotype_data) 

motif_and_phenotypes_count <- transpose_count_triplet_motifs_results(motif_and_phenotypes, 
																	 "type_ac", "type_bc", "total_count", "motif_type") 	

write.table(motif_and_phenotypes_count , file =  file.path( results_directory, observed_counts_phenotypes_file), 
			row.names = FALSE)

#########################################################

### Randomized triplet motifs analysis 10 times, example
### This uses the mclappy function to spead the calculations onto different cores. It will only work on Linux machines and not Windows.

edited_fitered_costanzo_stringent <- filtered_costanzo_stringent

colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "query_oln_id_edited"] <- "oln_id_a"
colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "array_oln_id_edited"] <- "oln_id_b"

### Run one trial to see whether it works or not 
# result_one_trial <- run_one_randomized_trial_compare_with_dataset(1, edited_fitered_costanzo_stringent, kinase_network_subset,
# 															      sbi_interactome, tf_network, num_iterations=100,
# 																  FUNCT=count_triplet_motifs_phenotypes,
# 																  min_num_neighbours=10, max_num_tries=10, mode="dd", use_rewired_interactions=FALSE,
# 																  is_rewire_genetic_interactions=TRUE,
# 																  filtered_phenotype_data)
# 			

# Need to have mc.set.seed = TRUE and mc.preschedule=FALSE as per user manual on ?mcparallel in the parallel library. See section on Random Numbers.
## "The behaviour with mc.set.seed = TRUE is different only if RNGkind("L'Ecuyer-CMRG") has been selected. Then each time a child is forked it is given the next stream (see nextRNGStream). So if you select that generator, set a seed and call mc.reset.stream just before the first use of mcparallel the results of simulations will be reproducible provided the same tasks are given to the first, second, ... forked process." -- user manual on ?mcparallel in the parallel library, section on Random Numbers

mc.reset.stream() 
randomized_counts_phenotypes_list <- mclapply ( X=1:number_of_randomized_trials, FUN= run_one_randomized_trial_compare_with_dataset,   
														  edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
														  num_iterations=num_iteration_rewire_network, 
														  count_triplet_motifs_phenotypes, 
														  min_num_neighbours=10, max_num_tries=10, mode="dd", use_rewired_interactions=FALSE,
											  			  is_rewire_genetic_interactions=TRUE,
														  filtered_phenotype_data,
														  mc.set.seed = TRUE , 
														  mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )

randomized_counts_phenotypes_table <- concat_motif_counts_list_into_table( randomized_counts_phenotypes_list) 

write.table ( randomized_counts_phenotypes_table, file = file.path( results_directory, randomized_counts_phenotypes_file), 
			  row.names = FALSE)

#########################################################

full_results_phenotypes <- get_full_results_table(motif_and_phenotypes_count, 
																randomized_counts_phenotypes_table, p.value=0.05) 

write.table ( full_results_phenotypes, file = file.path( results_directory, full_results_phenotypes_file), 
			  row.names = TRUE) 

#########################################################

