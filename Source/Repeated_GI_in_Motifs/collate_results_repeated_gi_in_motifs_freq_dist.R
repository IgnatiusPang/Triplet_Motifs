library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(purrr)

#####################################################################								 

### Graphic width and heigth for the boxplots of the randomization results. The observed value is shown in a dot. 
### This setting is for the poster presentation.


source( "./Common/parameters_file.R")

source( file.path(source_directory_common, "count_triplet_motifs_helper.R" ) )
source( file.path(source_directory_common, "concatenate_result_files.R") )
source( file.path(source_directory, 'Figures/paper_figures_helper.R') )

results_directory       <- file.path( base_directory, "Results/Bootstrap_p_values/Repeated_GI_in_Motifs_Freq_Dist" )
final_results_directory <- file.path ( results_directory, "Final_Results" )

if (! file.exists(final_results_directory)){
	  dir.create(final_results_directory)
}

################################################################################################################################################


convert_panel_group_name <- function( input) {
	
	if( input =="protein_complexes" | input ==  "Protein Complexes") {
		return ( "Protein Complexes")
	} else if ( input ==  "regulatory_triplets" | input==  "Regulatory Triplets" ){
		return( "Regulatory Triplets")
	} else if ( input == "signaling_triplets" | input == "Signaling Triplets" ) {
		return ( "Signaling Triplets")
	} else if ( input == "others" | input == "Others" ) {
		return ( "Others")
	}
	
}



################################################################################################################################################
# Counts of Triplet Motifs Costanzo et al. 2010 dataset

triplet_motifs_full_results_file     <- "full_results_counts_unique_gi_pairs_in_motifs.tab"
count_triplet_motifs_randomized_file <- "randomized_counts_unique_gi_pairs_in_motifs.tab" 

triplet_motifs_file_pattern      	 <- "^randomized_counts_unique_gi_pairs_in_motifs.tab"
observed_file_pattern        		 <- "Job_20/observed_counts_unique_gi_pairs_in_motifs.tab"
final_results_directory_pattern      <- "Final_Results"
number_of_randomized_samples         <- 68544 
p_values 							 <- 0.05

## Randomized data
randomized_table <- collate_randomization_result_files_into_one_table( results_directory, triplet_motifs_file_pattern,  
																	   final_results_directory_pattern, 
																  	   number_of_randomized_samples, id =TRUE )


randomized_table <- randomized_table %>%
	dplyr::mutate( motif_type = sapply( motif_type, convert_panel_group_name ) )


# collated_table <- collate_result_files_helper ( results_directory, triplet_motifs_full_results_file, count_triplet_motifs_randomized_file,
# 							  triplet_motifs_file_pattern, observed_file_pattern, final_results_directory_pattern,
# 							  number_of_randomized_samples, p_values)

#   counts file_number group motif_type total_count
# 1      2           1     1     others          60
# 2      3           1     1     others          18
# 3      4           1     1     others           3
# 4      5           1     1     others           4
# 5      6           1     1     others           3
# 6      7           1     1     others           4

## Observed Data
observed_data <- read.table ( file.path( results_directory, observed_file_pattern  ), header=TRUE)

observed_data <- observed_data %>%
	dplyr::mutate( motif_type = sapply( motif_type, convert_panel_group_name ) )

observed_data <-  dplyr::filter(observed_data, motif_type!="Others")



## Find the minimum and maximum number of overlapping triplets. 
min_max_counts_filter <-  observed_data %>% 
	dplyr::group_by(motif_type) %>%
	dplyr::summarise( minimum=min(counts), maximum=max(counts)) %>%
	dplyr::select( one_of( c( "motif_type", "minimum", "maximum"))) %>%
	dplyr::mutate( motif_type = sapply( motif_type, convert_panel_group_name ) )

## Get all possible combinations of file_number,  motif_type, group, counts
combinations <- purrr::cross_df(list(file_number=unique(randomized_table[, 'file_number'] ), 
				 group=unique(randomized_table[, 'group']),
				 motif_type=unique(randomized_table[, 'motif_type']), 
				 counts = 2:88)) %>%
				 dplyr::mutate( motif_type = sapply( motif_type, convert_panel_group_name ) )%>%
				 dplyr::left_join ( min_max_counts_filter, by="motif_type" ) %>%
				 dplyr::filter ( counts >= minimum & counts <= maximum) %>%
				 dplyr::arrange( file_number,  motif_type, group, counts) 

combinations <- combinations %>% dplyr::filter ( counts >= minimum & counts <= maximum)

### Use all the combinations possible to fill in counts that were observed with zero 
randomized_table <- combinations %>% 
	dplyr::left_join ( randomized_table ) %>%
	dplyr::mutate( total_count = ifelse( is.na(total_count), 0, total_count)) %>%
	dplyr::select ( -minimum, -maximum) %>%
	dplyr::mutate( motif_type = sapply( motif_type, convert_panel_group_name ) )

# Calculate adjusted p-values
statistical_analysis <- observed_data %>% 
	dplyr::mutate( motif_type = sapply( motif_type, convert_panel_group_name ) ) %>%
	dplyr::rename( observed_count = total_count) %>%
	dplyr::right_join( randomized_table,  by=c( "motif_type" = "motif_type", "counts" = "counts")) %>%
	dplyr::filter( motif_type != 'Others') %>%
	dplyr::mutate( observed_count =  ifelse( is.na(observed_count), 0, observed_count)) %>%
	dplyr::mutate( is_above = ifelse ( observed_count > total_count, 1, 0 ) ) %>%
	dplyr::group_by( motif_type, counts, observed_count) %>%
	dplyr::summarise( p_value = 1-  sum(is_above)/n(), sum_if_above = sum(is_above), num_entries = n()  ) %>%  ## Calculate raw p-values
	dplyr::select ( one_of ( c( "motif_type", "counts", "observed_count", "p_value") )) %>%
	dplyr::left_join ( min_max_counts_filter) %>%                          
	dplyr::mutate( adj_p_value = p_value * ( maximum - minimum + 1)) %>% ## Adjusted p-values
	dplyr::mutate( adj_p_value = ifelse ( adj_p_value > 1, 1, adj_p_value )   )  %>%
	dplyr::mutate( is_significant = ifelse ( adj_p_value < 0.05, TRUE, FALSE )) %>%
	dplyr::select ( one_of ( c( "motif_type", "counts", "observed_count", "p_value", "adj_p_value", "is_significant") )) %>%
	dplyr::rename(total_count = observed_count) %>%
	dplyr::mutate( is_significant_color = factor ( boolean_to_colour(is_significant), levels=c("red", "blue") )) %>%
	dplyr::mutate( total_count = total_count + 1) # %>%  ### Add 1 as a pseudo count for drawing in log scale
	# dplyr::filter( counts <= max_num_overlapping_triplets_to_show )  
     ## For the protein complexes class, there is a genetic interaction with up to ~80 overlapping triplets, so we want to trim the x-axis smaller.
	

write.table( randomized_table, file.path( results_directory, "Final_Results/overlapping_triplets_randomized_frequency_distribution.txt"))

write.table( statistical_analysis, file.path( results_directory, "Final_Results/overlapping_triplets_frequency_distribution.txt"))

#################################################################

