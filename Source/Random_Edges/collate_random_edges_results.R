library(tidyr)
library(dplyr)
library(ggplot2)
library(lazyeval)


options <- commandArgs(trailingOnly = TRUE)
source( "./Common/parameters_file.R")


# figure_file_suffix <- ".tiff"
figure_file_suffix <- ".pdf"

source( file.path ( source_directory, 'Random_Edges/collate_random_edges_results_helper.R')) 
source( file.path( source_directory, 'Figures/paper_figures_helper.R') )

### Directories
results_directory <- file.path ( results_directory, "Bootstrap_p_values/Random_Edges" ) 
final_results_directory <- file.path ( results_directory, "Final_Results")

create_directory_if_not_exists(results_directory)
create_directory_if_not_exists(final_results_directory)


### Size of resulting plot 
graphic_width  <- 10 
graphic_height <- 10
	
### Parameter space
PERCENTAGE_OF_ORIGINAL_NETWORK <- seq( 0.2, 1.3, 0.1 )
INNER_LOOP <- seq(1, 20, 1)

JOB_NUMBERS <- 0:239  ## Need to be upldated once the corrected results has been calculated by Katana
adjust_job_number <- 1 ## Add this number to the job number to access the correct position in the  PARAMS array

## Use purple background and create a subset of the triplet motifs for the poster
print_poster_results_purple_background <- FALSE

plot_background_colour <- "white"  
axis_text_colour       <- "black" 

if ( print_poster_results_purple_background== TRUE) {
	plot_background_colour <-  "white" # "#7E4E99" 
	axis_text_colour       <-  "black" # "#FFC000" 
}

##########################################################################################################

PARAMS <- c()

for( PERCENTAGE in PERCENTAGE_OF_ORIGINAL_NETWORK ) {
	
	for ( LOOP_COUNT in INNER_LOOP ) { 

		PARAMS <- c( PARAMS, PERCENTAGE) 
	}
}

# random_edges_after_randomization_counts_0.2_1_1_1_.tab
# random_edges_before_randomization_counts_0.2_1_1_1_.tab

prev_step_network_size <- 0 

before_randomization_file_list <- c()
after_randomization_file_list <- c() 

for ( job_number in JOB_NUMBERS) {

	proprotion_of_original_network_size <- PARAMS[job_number+adjust_job_number]

	
	before_randomization_file_list <- c(before_randomization_file_list, 
										paste ( "Job_", job_number, 
										 "/random_edges_before_randomization_counts_", 
										 proprotion_of_original_network_size, 
										 "_1_1_1_.tab", sep="") ) 
	
	after_randomization_file_list <- c( after_randomization_file_list, 
										paste ( "Job_", job_number, 
										"/random_edges_after_randomization_counts_", 
										proprotion_of_original_network_size, 
										"_1_1_1_.tab", sep="") ) 
}

before_randomization_collated_counts <- collate_different_datasets(  results_directory, before_randomization_file_list, 
																	 strsplit_pattern = '_', position_of_parameter = 7, 'Before')  
		
after_randomization_collated_counts <- collate_different_datasets( results_directory, after_randomization_file_list, 
																	strsplit_pattern = '_', position_of_parameter = 7, 'After')  
				
randomization_collated_counts <- rbind( before_randomization_collated_counts, after_randomization_collated_counts)

## Fix the naming of motif types 

randomization_collated_counts[, 'motif_type'] <- convert_triplet_motifs_name_to_paper_style( randomization_collated_counts[, 'motif_type'])



#####################################################################################################################################

stat_test_result_random_edges <- random_edges_statstical_testing(randomization_collated_counts, 
																 counts_percentage_threshold = 0.02, p_value_threshold= 0.05, 
																 num_of_motifs=15, num_random_trials=2000, 
																 experiment_label_observed='Before',
																 experiment_label_random='After' ) 

stat_test_result_random_edges[, 'motif_type'] <- convert_triplet_motifs_name_to_paper_style( stat_test_result_random_edges[, 'motif_type'])

write.table ( stat_test_result_random_edges, file=file.path ( final_results_directory, "stat_test_result_random_edges.txt"))

### The trend is that the number of triplet motifs will get less, but there are still likely to be significant differences with a random network. 

### Ignore the pairedness of the analyses

#####################################################################################################################################

# The motif types sorted by largest to lowest count
motif_types_sorted_by_counts_desc <- dplyr::filter ( randomization_collated_counts, parameter==1 & experiment=='Before' & motif_type != 'total' ) %>%
									 dplyr::select( one_of(c('experiment', 'motif_type', 'counts'))) %>% 
									 dplyr::arrange(desc(counts)) %>%
									 dplyr::distinct() %>% 
									 dplyr::select (one_of(c('motif_type'))) %>% 
									 as.data.frame() %>%
									 t() %>%
									 as.vector()

## Temporary code to filter data for Poster
if ( print_poster_results_purple_background == TRUE) {
	randomization_collated_counts <-dplyr::filter (randomization_collated_counts,  
												   motif_type %in% c( 'PP', 'TUTU', 'PKU', 'KUKU', 'TDP', 'PKD')) 
	
	motif_types_sorted_by_counts_desc <- motif_types_sorted_by_counts_desc[motif_types_sorted_by_counts_desc %in% c( 'PP', 'TUTU', 'PKU', 
																													 'KUKU', 'TDP', 'PKD')]	
	stat_test_result_random_edges <- dplyr::filter (stat_test_result_random_edges,  motif_type %in% c( 'PP', 'TUTU', 'PKU', 
																									   'KUKU', 'TDP', 'PKD')) 
}

### Get statistical significance of results
background_shading_table <- dplyr::select(stat_test_result_random_edges, one_of(c('parameter', 'motif_type', 'is_significant')))
background_shading_table$experiment <-1
background_shading_table$counts <-1
background_shading_table[, 'is_significant'] <- as.factor(background_shading_table[, 'is_significant'] )

## Change the order in which the motif types are ordered when they appear as rows in the facet grid
background_shading_table[, 'motif_type'] <- factor (  background_shading_table[, 'motif_type'], levels = motif_types_sorted_by_counts_desc ) 
randomization_collated_counts[, 'motif_type'] <- factor ( randomization_collated_counts[,'motif_type'], levels = motif_types_sorted_by_counts_desc ) 


#####################################################################################################################################

write.table ( randomization_collated_counts, file=file.path ( final_results_directory, "randomization_collated_counts.txt"))

#####################################################################################################################################

