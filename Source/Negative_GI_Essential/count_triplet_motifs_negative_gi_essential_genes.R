### Script: count_triplet_motifs_negative_gi_essential_genes.R
### Author: Ignatius Pang 
### Date: 27-5-2016
### Description: Count the observed number of essential genes in triplet motifs. 
### Randomize the respective networks and count the number of esential genes in triplet motifs. Keep the essential genes the same, only randomize the edges.
### Repeat randomization 2000 times.
### Obtain the bootstrap p-value

#########################################################
# Source location

# psql sbi_triplet_motifs
# cd '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/'
# 

#########################################################
### Parameters 
options <- commandArgs(trailingOnly = TRUE)


### Source the parameters file here 

if ( length(options) != 0 
	 & ( options[1] == 'katana' | options[1] == 'clive')  ) {
	source( "./parameters_file.R" )

} else {
	source( "/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Common/parameters_file.R")
}


#### Local Parameters
if (is_run_locally) {
	results_directory <- paste ( base_directory, "Results/Bootstrap_p_values_temp/Negative_GI_Essential/",  sep="")
}

observed_a_and_b_essential_gene_counts 	<- "observed_a_and_b_essential_gene_counts.tab"
observed_c_essential_gene_counts 	<- "observed_c_essential_gene_counts.tab"

randomized_count_a_and_b_table_file 	<- "randomized_results_table_count_a_and_b.tab"
randomized_count_c_table_file 		<- "randomized_results_table_count_c.tab"

output_full_results_table_count_a_and_b	<- "full_results_table_count_a_and_b.tab"		  
output_full_results_table_count_c  	<- "full_results_table_count_c.tab"


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


## Insert code for counting essential genes here 
## A table contiaining the counts of the number of essential genes within each type of triplet motifs
count_essential_genes_in_motifs <- count_triplet_motifs_against_gene_list( triplet_motifs_costanzo, essential_genes)


results_count_a_and_b <- transpose_count_triplet_motifs_results(count_essential_genes_in_motifs[,c("type_ac", "type_bc", "count_a_and_b")], 
												 "type_ac", "type_bc", "count_a_and_b", "motif_type") 


results_count_c <- transpose_count_triplet_motifs_results(count_essential_genes_in_motifs[,c("type_ac", 
												"type_bc", "count_c")], 
												"type_ac", "type_bc", "count_c", "motif_type") 

write.table(results_count_a_and_b , file =  paste( results_directory, observed_a_and_b_essential_gene_counts, sep=""), 
			row.names = FALSE)

write.table(results_count_c , file =  paste( results_directory, observed_c_essential_gene_counts, sep=""), 
			row.names = FALSE)


#########################################################

### Randomized triplet motifs analysis 10 times, example
### This uses the mclappy function to spead the calculations onto different cores. It will only work on Linux machines and not Windows.

edited_fitered_costanzo_stringent <- filtered_costanzo_stringent

colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "query_oln_id_edited"] <- "oln_id_a"
colnames( edited_fitered_costanzo_stringent)[colnames( edited_fitered_costanzo_stringent) == "array_oln_id_edited"] <- "oln_id_b"

### Run one trial to see whether it works or not 
# result_one_trial <- run_one_randomized_trial_compare_with_dataset(1, edited_fitered_costanzo_stringent, kinase_network_subset,
																  # sbi_interactome, tf_network,
																  # num_iterations=10, count_triplet_motifs_against_gene_list, 
																  # min_num_neighbours=10, max_num_tries=10, mode="dd",
																  # essential_genes)


# Need to have mc.set.seed = TRUE and mc.preschedule=FALSE as per user manual on ?mcparallel in the parallel library. See section on Random Numbers.
## "The behaviour with mc.set.seed = TRUE is different only if RNGkind("L'Ecuyer-CMRG") has been selected. Then each time a child is forked it is given the next stream (see nextRNGStream). So if you select that generator, set a seed and call mc.reset.stream just before the first use of mcparallel the results of simulations will be reproducible provided the same tasks are given to the first, second, ... forked process." -- user manual on ?mcparallel in the parallel library, section on Random Numbers

mc.reset.stream() 
list_of_randomized_triplet_motif_counts <- mclapply ( X=1:number_of_randomized_trials, FUN= run_one_randomized_trial_compare_with_dataset,   
													  edited_fitered_costanzo_stringent, kinase_network_subset, sbi_interactome, tf_network, 
													  num_iterations=num_iteration_rewire_network, 
													  count_triplet_motifs_against_gene_list, 
													  min_num_neighbours=10, max_num_tries=10, mode="dd",
													  essential_genes,
													  mc.set.seed = TRUE , 
													  mc.preschedule=FALSE, mc.cores=number_of_cores_to_use )

													  
### Count the number of essential genes among gene A and gene B of triplet motifs. Repeated for many randomized networks.
count_a_and_b_list <-  lapply ( list_of_randomized_triplet_motif_counts, function(x) {return( x[,c("type_ac", "type_bc", "count_a_and_b")] )} )
													  
randomized_count_a_and_b_table <- concat_motif_counts_list_into_table( count_a_and_b_list) 

write.table ( randomized_count_a_and_b_table, file = paste( results_directory, randomized_count_a_and_b_table_file, sep=""), 
			  row.names = FALSE)

			  
### Count the number of essential genes among gene A and gene B of triplet motifs. Repeated for many randomized networks.
count_c_list <-   lapply ( list_of_randomized_triplet_motif_counts, function(x) {return( x[,c("type_ac", "type_bc", "count_c")] )} ) 													
													  
randomized_count_c_table <- concat_motif_counts_list_into_table( count_c_list) 
			  
write.table ( randomized_count_c_table, file = paste( results_directory, randomized_count_c_table_file, sep=""), 
			  row.names = FALSE)			  

#########################################################


# Write the output tables, use the p-value for two-sided test for 'count gene A and gene B' 
final_results_table_count_a_and_b <- get_full_results_table(results_count_a_and_b, randomized_count_a_and_b_table, p.value=0.05) 

write.table ( final_results_table_count_a_and_b, file = paste( results_directory, output_full_results_table_count_a_and_b, sep=""), 
			  row.names = TRUE) 

# Write the output tables, use the p-value for two-sided test for 'count gene C' 
final_results_table_count_c <- get_full_results_table(results_count_c, randomized_count_c_table, p.value=0.05) 

write.table ( final_results_table_count_c, file = paste( results_directory, output_full_results_table_count_c, sep=""), 
			  row.names = TRUE) 



