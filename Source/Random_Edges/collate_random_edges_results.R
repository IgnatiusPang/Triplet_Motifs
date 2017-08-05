library(tidyr)
library(dplyr)
library(ggplot2)
library(lazyeval)


options <- commandArgs(trailingOnly = TRUE)
source( "./Common/parameters_file.R")


source( file.path ( source_directory, 'Random_Edges/collate_random_edges_results_helper.R')) 
source( file.path( source_directory, 'Figures/reshape_triplet_motifs_helper.R') )

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
print_poster_results_purple_background <- TRUE

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

write.table ( randomization_collated_counts, file=file.path ( final_results_directory, "randomization_collated_counts.txt"))


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
motif_types_sorted_by_counts_desc <- filter ( randomization_collated_counts, parameter==1 & experiment=='Before' & motif_type != 'total' ) %>%
									select( one_of(c('experiment', 'motif_type', 'counts'))) %>% 
									arrange(desc(counts)) %>%
									distinct() %>% 
										select (one_of(c('motif_type'))) %>% 
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
temp <- dplyr::select(stat_test_result_random_edges, one_of(c('parameter', 'motif_type', 'is_significant')))
temp$experiment <-1
temp$counts <-1
temp[, 'is_significant'] <- as.factor(temp[, 'is_significant'] )

## Change the order in which the motif types are ordered when they appear as rows in the facet grid
temp[, 'motif_type'] <- factor (  temp[, 'motif_type'], levels = motif_types_sorted_by_counts_desc ) 
randomization_collated_counts[, 'motif_type'] <- factor ( randomization_collated_counts[,'motif_type'], levels = motif_types_sorted_by_counts_desc ) 

# removal or addition of nodes to network
randomization_collated_counts %>%
	dplyr::filter (  motif_type != 'total' ) %>%
	ggplot (  aes( experiment, counts) ) +
	geom_boxplot()  +
	geom_rect(data = temp ,aes(  fill=is_significant ),
			  xmin = -Inf,xmax = Inf,
			  ymin = -Inf,ymax = Inf,alpha = 0.3) +
	scale_fill_manual( values=c('blue', 'red')) +
	facet_grid(   motif_type ~ parameter , scales="free") + 
	labs( title = "Removal or Addition of Edges to Network") +
	xlab( "Before or After Network Randomization") +
	ylab ('Counts') +
	theme( plot.background = element_rect(fill = plot_background_colour, color=plot_background_colour), 
		   axis.text = element_text(colour = axis_text_colour),
		   axis.title = element_text(colour = axis_text_colour),
			strip.text.y = element_text(face = "italic"),
			  plot.title = element_text(hjust = 0.5, size = rel(1), colour = axis_text_colour),
			  	legend.position="none")


if ( print_poster_results_purple_background == FALSE) {
	
ggsave(file.path(final_results_directory, "randomization_collated_results.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )

ggsave(file.path(final_results_directory, "randomization_collated_results_15x15.tiff"), plot=last_plot(), width=15, height=15 )

}

if ( print_poster_results_purple_background == TRUE) {
	
	final_results_directory <-	"/media/z3371724/PostDoc/2016/Triplet_Motifs/Results/Poster/Figures_Purple_Background/"
	
	ggsave(file.path(final_results_directory, "randomization_collated_results.tiff"), plot=last_plot(), width=10, height=5 )

}


#####################################################################################################################################


