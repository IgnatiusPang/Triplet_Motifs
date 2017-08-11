### Script: count_triplet_motifs_random_edges.R
### Author: Ignatius Pang 
### Date: 20-5-2016
### Description:Randomly add or remove edges to change the size of the network based on the original network size, count the number of motifs. Randomize the network and count the number triplet motifs. Repeat the above steps 2000 times.
### Obtain the bootstrap p-value

library(igraph)
library(dplyr)
library(lazyeval)
library(parallel)

sessionInfo()

#########################################################
# Source location

# psql sbi_triplet_motifs
# Rscript --vanilla count_triplet_motifs_random_edges.R local 0.5 1 1 1  

#########################################################
### Parameters 
options <- commandArgs(trailingOnly = TRUE)

### Source the parameters file here 
if ( length(options) != 0  )  { 
	if ( options[1] == 'katana' | options[1] == 'clive')   {
		source( "./parameters_file.R" )
	} else {
		print ('Load parameter file locally.')
		source( "./Common/parameters_file.R")
	}
	
} else {
	print ('Load parameter file locally.')
	source( "./Common/parameters_file.R")
}

source_directory_random_edges <- "./"

if (is_run_locally) {
	results_directory <- file.path( results_directory, "Bootstrap_p_values_temp/Random_Edges")
	source_directory_random_edges <- file.path( source_directory, "Random_Edges")
	
	create_directory_if_not_exists(results_directory)
}

source( file.path (source_directory_random_edges, "random_edges_helper.R" ))

### Local Parameters

### Randomly add or remove edges to change the size of the network based on the original network size
proportion_of_original_network_size_gi_network      <- 0.5 # Default no change at all
proportion_of_original_network_size_kinase_network  <- 1
proportion_of_original_network_size_ppi_network     <- 1 
proportion_of_original_network_size_tf_network      <- 1


if ( length(options) != 0 
	 & ( options[1] == 'katana' | options[1] == 'clive') ) {
	
	if ( length( options ) != 7 ) {

		 stop_program_notification_string <- paste('Not enough command line arguments.', 
		 									 'Rscript --vanilla <name of script.R> <katana | clive | local> <array job id>' , 
		 									 '<Number of Random Trials> <gi network proportion>',
											 '<kinase network proportion> <ppi network proportion> <tf network proportion>'    )
		 
		 stop( stop_program_notification_string)
		 
	} else {
		
		 proportion_of_original_network_size_gi_network     <- as.numeric (options[4] )
		 proportion_of_original_network_size_kinase_network <- as.numeric (options[5] )
		 proportion_of_original_network_size_ppi_network    <- as.numeric (options[6] )
		 proportion_of_original_network_size_tf_network     <- as.numeric (options[7] )
		 
	}
}

print ( paste ( "proportion_of_original_network_size_gi_network =", proportion_of_original_network_size_gi_network))

before_randomization_counts_file <- paste( "random_edges_before_randomization_counts_", proportion_of_original_network_size_gi_network, "_", 
										   proportion_of_original_network_size_kinase_network, "_", 
										   proportion_of_original_network_size_ppi_network, "_", 
										   proportion_of_original_network_size_tf_network, "_", 
										    ".tab", sep="")
										   					
after_randomization_counts_file  <- paste( "random_edges_after_randomization_counts_", proportion_of_original_network_size_gi_network, "_", 
										   proportion_of_original_network_size_kinase_network, "_", 
										   proportion_of_original_network_size_ppi_network, "_",  
										   proportion_of_original_network_size_tf_network, "_", ".tab", sep="")


#########################################################

run_add_or_remove_edges_once <- function(iteration_holder, filtered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
										 proportion_of_original_network_size_gi_network = 1,
										 proportion_of_original_network_size_kinase_network = 1,
										 proportion_of_original_network_size_ppi_network = 1,
										 proportion_of_original_network_size_tf_network = 1,
										 num_iterations=NULL) {
  	
	##
	filtered_costanzo_stringent <- randomly_add_or_remove_edges_from_network (  filtered_costanzo_stringent,  proportion_of_original_network_size_gi_network,
																	"query_oln_id_edited", "array_oln_id_edited",
																	directed = FALSE, 											 			
																	simple_network = FALSE,
																	maximum_num_of_tries = 10 )
	

	temp_vector <- rep ( 0, count(filtered_costanzo_stringent)[[1]])

	filtered_costanzo_stringent <- cbind ( filtered_costanzo_stringent, genetic_interaction_score= temp_vector, p_value=temp_vector, std_dev=temp_vector)

    ##
	tf_network <- randomly_add_or_remove_edges_from_network (  tf_network,  proportion_of_original_network_size_tf_network,
                        												 "regulator_oln_id", "target_oln_id",
                        												 directed = TRUE,
                        												 maximum_num_of_tries = 10 )
	

	##	
	kinase_network_subset <- randomly_add_or_remove_edges_from_network (  kinase_network_subset,  proportion_of_original_network_size_kinase_network,
	                                                          "kinase_oln_id", "target_oln_id",
	                                                          directed = TRUE,
	                                                          maximum_num_of_tries = 10 )
	
	
	##
	sbi_interactome <- randomly_add_or_remove_edges_from_network (  sbi_interactome,  proportion_of_original_network_size_ppi_network, 
														"oln_id_a", "oln_id_b",
														directed = FALSE, 
														maximum_num_of_tries = 10 )
	

#########################################################

### Createa a table that include all kinase-substrate, protein-protein and transcription factor-gene interactions
### Do not contain genetic interactions

tf_network_collated <- collate_interactions_from_both_direction_tbl_df(tbl_df(tf_network), "regulator_oln_id", "target_oln_id",
                                                                       "transcription factor-target down", "transcription factor-target up",
                                                                       "td", "tu",
                                                                       "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev")


kinase_network_collated <- 	tbl_df(kinase_network_subset) %>%
  collate_interactions_from_both_direction_tbl_df( "kinase_oln_id", "target_oln_id",
                                                   "kinase-substrate down", "kinase-substrate up",
                                                   "kd", "ku",
                                                   "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev")

sbi_interactome_collated <- collate_interactions_from_both_direction_tbl_df( tbl_df(sbi_interactome), "oln_id_a", "oln_id_b",
                                                                             "protein-protein", "protein-protein",
                                                                             "p", "p",
                                                                             "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev")


# ### All interactions except for genetic interactions
interactions_combined <- dplyr::bind_rows(tf_network_collated, kinase_network_collated, sbi_interactome_collated)


#########################################################

#### Get the number of triplet motifs before randomization

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

# results_count_triplet_motifs_costanzo <- transpose_count_triplet_motifs_results  (results_count_triplet_motifs_costanzo,  "type_ac", "type_bc", "my_total_count", "motif_type") 

# print ( results_count_triplet_motifs_costanzo)

#########################################################

#### Get the number of triplet motifs after randomization

### Randomized each network,
### This uses the mclappy function to spead the calculations onto different cores. It will only work on Linux machines and not Windows.

edited_fitered_costanzo_stringent <- filtered_costanzo_stringent

colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "query_oln_id_edited"] <- "oln_id_a"
colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "array_oln_id_edited"] <- "oln_id_b"


### Run one trial to see whether it works or not 
result_one_trial <- run_one_randomized_trial_compare_with_dataset(1, edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network,
 																   num_iterations=num_iteration_rewire_network, count_triplet_motifs)


# result_one_trial <- transpose_count_triplet_motifs_results  (result_one_trial,  "type_ac", "type_bc", "total_count", "motif_type")


# print( result_one_trial)

return( list( before_randomization= results_count_triplet_motifs_costanzo, after_randomization=result_one_trial) ) 

}

# 
# test_one_result <- run_add_or_remove_edges_once (1, filtered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
# 												 proportion_of_original_network_size_gi_network = proportion_of_original_network_size_gi_network ,
# 												 proportion_of_original_network_size_kinase_network = proportion_of_original_network_size_kinase_network ,
# 												 proportion_of_original_network_size_ppi_network = proportion_of_original_network_size_ppi_network ,
# 												 proportion_of_original_network_size_tf_network = proportion_of_original_network_size_tf_network ,
# 												 num_iterations = num_iteration_rewire_network)
# 
# print (test_one_result )
# 	

mc.reset.stream()
list_of_randomized_triplet_motif_counts <- mclapply ( X=1:number_of_randomized_trials, FUN= run_add_or_remove_edges_once,
						filtered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network,
						proportion_of_original_network_size_gi_network = proportion_of_original_network_size_gi_network ,
						proportion_of_original_network_size_kinase_network = proportion_of_original_network_size_kinase_network ,
						proportion_of_original_network_size_ppi_network = proportion_of_original_network_size_ppi_network ,
						proportion_of_original_network_size_tf_network = proportion_of_original_network_size_tf_network ,
						num_iterations = num_iteration_rewire_network,
						mc.set.seed = TRUE ,
						mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )

if ( 'before_randomization' %in% names(list_of_randomized_triplet_motif_counts[[1]] )  ) {
	
	### Count the number of periodic genes among gene A and gene B of triplet motifs. Repeated for many randomized networks.
	before_randomization_list <-  lapply ( list_of_randomized_triplet_motif_counts, function(x) {return( x$before_randomization )} )
	
	before_randomization_table <- concat_motif_counts_list_into_table( before_randomization_list)
	
	write.table ( before_randomization_table, file = file.path( results_directory, before_randomization_counts_file),
				  row.names = FALSE)

	after_randomization_list <-  lapply ( list_of_randomized_triplet_motif_counts, function(x) {return( x$after_randomization )} )
	
	after_randomization_table <- concat_motif_counts_list_into_table( after_randomization_list)
	
	write.table ( after_randomization_table, file = file.path( results_directory, after_randomization_counts_file),
				  row.names = FALSE)
	
} else {
	
	print ( list_of_randomized_triplet_motif_counts)
	
}

#########################################################

