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

	results_directory <- "/media/z3371724/PostDoc/2016/Triplet_Motifs/Results/Bootstrap_p_values/Protein_Complexes/"
}

observed_counts_babu_protein_complex_file        <- "observed_counts_babu_protein_complex_in_tri_motifs.tab"
observed_counts_benschop_protein_complex_file    <- "observed_counts_benschop_protein_complex_in_tri_motifs.tab"

randomized_counts_babu_protein_complex_file 	 <- "randomized_counts_babu_protein_complex_in_tri_motifs.tab"
randomized_counts_benschop_protein_complex_file  <- "randomized_counts_benschop_protein_complex_in_tri_motifs.tab"
			  
full_results_babu_protein_complex_file 		     <- "full_results_babu_protein_complex_in_tri_motifs.tab"
full_results_benschop_protein_complex_file  	 <- "full_results_table_benschop_protein_complex_in_tri_motifs.tab"

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
#### Analysis on Babu membrane protein complexes dataset	
# babu_motif_and_protein_complexes <- count_triplet_motifs_in_protein_complexes ( triplet_motifs_costanzo, babu_protein_complexes) 
	

# babu_protein_complex_count <- transpose_count_triplet_motifs_results(babu_motif_and_protein_complexes, 
																# "type_ac", "type_bc", "total_count", "motif_type") 	

# write.table(babu_protein_complex_count , file =  paste( results_directory, observed_counts_babu_protein_complex_file, sep=""), 
			# row.names = FALSE)


##################################################################################################################################################
### Helper analysis to get statistics for use in paper 
# benschop_motif_and_protein_complexes_temp <- count_triplet_motifs_in_protein_complexes_helper ( triplet_motifs_costanzo, benschop_protein_complexes)
# 
# benschop_motif_and_protein_complexes_temp <- as.data.frame( benschop_motif_and_protein_complexes_temp)
# 
# # maps to how many protein complexes
# length( unique( benschop_motif_and_protein_complexes_temp[, "complex_id_a"]) )
# 
# ## how many proteins among motifs in protein complexes
# 
# length( unique( c( benschop_motif_and_protein_complexes_temp[, "oln_id_a"]
# 					, benschop_motif_and_protein_complexes_temp[, "oln_id_b"]
# 					, benschop_motif_and_protein_complexes_temp[, "oln_id_c"] ) )  )
# 
# # maps to how many protein complexes for pp motif
# 
# dplyr::filter ( benschop_motif_and_protein_complexes_temp, type_ac == 'p' & type_bc == 'p') %>%
# 	dplyr::select ( one_of ( c("complex_id_a"))) %>%
# 	dplyr::distinct() %>%
# 	dplyr::count()
# 
# ## how many proteins among pp motifs in protein complexes
# benschop_motif_pp <- dplyr::filter ( benschop_motif_and_protein_complexes_temp, type_ac == 'p' & type_bc == 'p')
# 
# ( dplyr::select (benschop_motif_pp,  one_of ( c("oln_id_a"))) %>% rename( oln_id = oln_id_a) )  %>%
# dplyr::union( dplyr::select (benschop_motif_pp,  one_of ( c("oln_id_b"))) %>% rename( oln_id = oln_id_b) ) %>%
# dplyr::union( dplyr::select (benschop_motif_pp,  one_of ( c("oln_id_c"))) %>% rename( oln_id = oln_id_c) )  %>%
# dplyr::distinct() %>%
# dplyr::count()	

##################################################################################################################################################
#### Analysis on Benschop core protein complexes dataset	


benschop_motif_and_protein_complexes <- count_triplet_motifs_in_protein_complexes ( triplet_motifs_costanzo, benschop_protein_complexes) 


benschop_protein_complex_count <- transpose_count_triplet_motifs_results(benschop_motif_and_protein_complexes, 
																	 "type_ac", "type_bc", "total_count", "motif_type") 	

write.table(benschop_protein_complex_count , file =  paste( results_directory, observed_counts_benschop_protein_complex_file, sep=""), 
			row.names = FALSE)


#########################################################

### Randomized triplet motifs analysis 10 times, example
### This uses the mclappy function to spead the calculations onto different cores. It will only work on Linux machines and not Windows.

edited_fitered_costanzo_stringent <- filtered_costanzo_stringent

colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "query_oln_id_edited"] <- "oln_id_a"
colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "array_oln_id_edited"] <- "oln_id_b"



### Run one trial to see whether it works or not 
# result_one_trial <- run_one_randomized_trial_protein_complexes(1, edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, babu_protein_complexes, 
# 			num_iterations=100)


# result_one_trial <- run_one_randomized_trial_compare_with_dataset(1, edited_fitered_costanzo_stringent, 
# 																  kinase_network_subset, sbi_interactome, tf_network, 
#																  num_iterations=100,	
# 																  count_triplet_motifs_in_protein_complexes, babu_protein_complexes )
# 																  


# Need to have mc.set.seed = TRUE and mc.preschedule=FALSE as per user manual on ?mcparallel in the parallel library. See section on Random Numbers.
## "The behaviour with mc.set.seed = TRUE is different only if RNGkind("L'Ecuyer-CMRG") has been selected. Then each time a child is forked it is given the next stream (see nextRNGStream). So if you select that generator, set a seed and call mc.reset.stream just before the first use of mcparallel the results of simulations will be reproducible provided the same tasks are given to the first, second, ... forked process." -- user manual on ?mcparallel in the parallel library, section on Random Numbers

mc.reset.stream() 
# randomized_counts_babu_protein_complex_list <- mclapply ( X=1:number_of_randomized_trials, FUN= run_one_randomized_trial_compare_with_dataset,   
													  # edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
													  # num_iterations=num_iteration_rewire_network, 
													  # count_triplet_motifs_in_protein_complexes, babu_protein_complexes,
													  # mc.set.seed = TRUE , 
													  # mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )

													  
											  
# randomized_counts_babu_protein_complex_table <- concat_motif_counts_list_into_table( randomized_counts_babu_protein_complex_list) 

# write.table ( randomized_counts_babu_protein_complex_table, file = paste( results_directory, randomized_counts_babu_protein_complex_file, sep=""), 
			  # row.names = FALSE)

			  

###
randomized_counts_benschop_protein_complex_list <- mclapply ( X=1:number_of_randomized_trials, FUN= run_one_randomized_trial_compare_with_dataset,   
														  edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
														  num_iterations=num_iteration_rewire_network, 
														  count_triplet_motifs_in_protein_complexes, 
														  min_num_neighbours=10, max_num_tries=10, mode="dd",
														  benschop_protein_complexes,
														  mc.set.seed = TRUE , 
														  mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )


randomized_counts_benschop_protein_complex_table <- concat_motif_counts_list_into_table( randomized_counts_benschop_protein_complex_list) 

write.table ( randomized_counts_benschop_protein_complex_table, file = paste( results_directory, randomized_counts_benschop_protein_complex_file, sep=""), 
			  row.names = FALSE)

#########################################################
# 
# temp_results_directory <- "/media/z3371724/PostDoc/2016/Triplet_Motifs/Results/Bootstrap_p_values/Protein_Complexes/"
# 
# benschop_protein_complex_count 					<- read.table ( paste( temp_results_directory, 
# 															 			 "observed_counts_benschop_protein_complex_in_tri_motifs.tab" , sep="" ) , header=TRUE)
# 
# randomized_counts_benschop_protein_complex_table  <- read.table ( paste( temp_results_directory, 
# 																		 "randomized_counts_benschop_protein_complex_in_tri_motifs.tab" , sep="" ) , header=TRUE)

# # Write the output tables, use the p-value for two-sided test for the Babu membrane protein complexes dataset
# full_results_babu_protein_complex <- get_full_results_table(babu_protein_complex_count, 
															# randomized_counts_babu_protein_complex_table, p.value=0.05) 

# write.table ( full_results_babu_protein_complex, file = paste( results_directory, full_results_babu_protein_complex_file, sep=""), 
			  # row.names = TRUE) 

# Write the output tables, use the p-value for two-sided test for the Benschop protein complexes dataset
full_results_benschop_protein_complex <- get_full_results_table(benschop_protein_complex_count, 
																randomized_counts_benschop_protein_complex_table, p.value=0.05) 

write.table ( full_results_benschop_protein_complex, file = paste( results_directory, full_results_benschop_protein_complex_file, sep=""), 
			  row.names = TRUE) 

#########################################################

