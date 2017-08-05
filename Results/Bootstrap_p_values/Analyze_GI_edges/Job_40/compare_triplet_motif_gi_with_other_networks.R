### Script: compare_triplet_motif_gi_with_other_networks.R
### Author: Ignatius Pang 
### Date: 22-8-2016
### Description: Compare the genetic interaction edge in triplet motifs with other types of interactions. 

#########################################################
# Source location

# psql sbi_triplet_motifs
# cd '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/'
# 

#########################################################
### Parameters 
### Global Parameters 

options <- commandArgs(trailingOnly = TRUE)

### Source the header file here 

if ( length(options) != 0 
	 & ( options[1] == 'katana' | options[1] == 'clive')  ) {
	source( "./parameters_file.R" )

} else {
	source( "/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Common/parameters_file.R")
}

### Local Parameters
observed_results_table_file   <- "observed_results_table_file.tab"
randomized_results_table_file <- "randomized_results_table.tab"
# output_enrichment_p_value     <- "enrichment_p_value_fixed.tab"
# output_depletion_p_value      <- "depletion_p_value_fixed.tab"
output_full_results_table     <- "full_results_table.tab"

if (is_run_locally) {
	
	results_directory <- "/media/z3371724/PostDoc/2016/Triplet_Motifs/Results/Bootstrap_p_values/Analyze_GI_edges/"
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

count_triplet_motifs_compare_with_other_networks <- function (triplet_motifs, type_ac, type_bc, interaction_type_abbrev, my_counts, total_count ) {
	
	a_ge_b <- filter_( triplet_motifs, paste(type_ac, ">=", type_bc) ) %>%
		group_by_( .dots=lapply ( c(type_ac, type_bc, interaction_type_abbrev), as.name)  ) %>%
		summarise_( .dots= list(counts=interp( ~sum(x), x=as.name(my_counts))) )
	
	columns_to_select <- c(type_ac, type_bc, interaction_type_abbrev, "counts")
	
	b_gt_a <- filter_( triplet_motifs, paste(type_bc, ">", type_ac) ) %>%
		group_by_( .dots=lapply ( c(type_ac, type_bc, interaction_type_abbrev), as.name) ) %>%
		summarise_( .dots= list(counts=interp( ~sum(x), x=as.name(my_counts))) ) %>%
		dplyr::rename_( .dots=setNames(list(type_ac, type_bc, interaction_type_abbrev), c(type_bc, type_ac, interaction_type_abbrev) ) ) %>%
		select( one_of( columns_to_select ))
	
	b_gt_a <- as.data.frame( b_gt_a )

	b_gt_a[,interaction_type_abbrev] <- sapply(b_gt_a[,interaction_type_abbrev],
												fix_edge_attribute_direction)

	b_gt_a <- dplyr::tbl_df(b_gt_a)
	
	triplet_motif_counts <- dplyr::union( a_ge_b, b_gt_a)  %>%
		group_by_( .dots=lapply ( c(type_ac, type_bc, interaction_type_abbrev), as.name) ) %>%	
		summarise( temp_counts=sum(counts))  %>%
		dplyr::rename_( .dots=setNames(list("temp_counts"), c(total_count) ) ) %>%
		arrange_( .dots=c(type_ac, type_bc) )
	
	return( triplet_motif_counts)
}

#########################################################


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

compare_gi_with_other_networks_final_results <- count_triplet_motifs_compare_with_other_networks_second_step( filtered_costanzo_stringent, interactions_combined, 
																											  "oln_id_a", "oln_id_b",
																											  "type_ac", "type_bc", "interaction_type_abbrev", 
																											  "count", "my_total_count")

write.table ( compare_gi_with_other_networks_final_results, file=paste(results_directory, 
																	   observed_results_table_file, 
																	   sep=""), row.names = FALSE)

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



### Concatenate the results together 
rbind_list_elements_recursively <- function ( x) {
	
	if (length(x)== 1 ) {
		return ( x[[1]])
	} else {
	
		return ( rbind ( x[[1]],  rbind_list_elements_recursively(x[2:length(x)]) ))
	}
}


randomized_results_table <- rbind_list_elements_recursively(list_of_randomized_triplet_motif_counts)

write.table ( randomized_results_table, file = paste( results_directory, randomized_results_table_file, sep=""), 
			  row.names = FALSE)

#########################################################


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
