library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

#####################################################################								 

### Graphic width and heigth for the boxplots of the randomization results. The observed value is shown in a dot. 
### This setting is for the poster presentation.

base_directory <- "/home/ignatius/PostDoc/2016/Triplet_Motifs/"
figures_results_directory <- "/home/ignatius/PostDoc/2016/Triplet_Motifs/Results/Poster/Figures"
	
source_directory <- paste( base_directory, "Source/", sep="")
source_directory_common <- paste( source_directory, "Common/",  sep="")
source( paste(source_directory_common, "count_triplet_motifs_helper.R", sep="") )
source(  paste(source_directory_common, "concatenate_result_files.R", sep="") )

results_directory                    <- paste( base_directory, "Results/Bootstrap_p_values/Repeated_GI_in_Motifs/", sep="")
final_results_directory              <- paste ( results_directory, "Final_Results/", sep="")

	if (! file.exists(final_results_directory)){
		dir.create(final_results_directory)

	}

################################################################################################################################################
# Counts of Triplet Motifs Costanzo et al. 2010 dataset

triplet_motifs_full_results_file     <- "full_results_counts_unique_gi_pairs_in_motifs.tab"
count_triplet_motifs_randomized_file <- "randomized_counts_unique_gi_pairs_in_motifs.tab" 

triplet_motifs_file_pattern      	 <- "^randomized_counts_unique_gi_pairs_in_motifs.tab"
observed_file_pattern        		 <- "Job_20/observed_counts_unique_gi_pairs_in_motifs.tab"
final_results_directory_pattern      <- "Final_Results"
number_of_randomized_samples         <- 2000
p_values 							 <- 0.05

collated_table <- collate_result_files_helper ( results_directory, triplet_motifs_full_results_file, count_triplet_motifs_randomized_file,
							  triplet_motifs_file_pattern, observed_file_pattern, final_results_directory_pattern,
							  number_of_randomized_samples, p_values)



################################################################################################

repeated_gi_full     <- read.table ( file.path( computational_results_directory,
											"Repeated_GI_in_Motifs/Final_Results/full_results_counts_unique_gi_pairs_in_motifs.tab"), header =TRUE)
repeated_gi_random   <- read.table ( file.path( computational_results_directory,
											"Repeated_GI_in_Motifs/Final_Results/randomized_counts_unique_gi_pairs_in_motifs.tab"), header =TRUE)
repeated_gi_observed <- read.table ( file.path( computational_results_directory,
											"Repeated_GI_in_Motifs/Job_1/observed_counts_unique_gi_pairs_in_motifs.tab"), header =TRUE)

observed_data <- repeated_gi_observed %>% tidyr::gather( "key", "value", 1:4) %>% dplyr::filter( key != 'others')
random_data <- repeated_gi_random %>% tidyr::gather( "key", "value", 1:4) %>% dplyr::filter( key != 'others')
p_values_data <- repeated_gi_full %>% dplyr::filter( motif_type != 'others')
	
# print_box_plot_observed_vs_random(repeated_gi_observed, repeated_gi_random, repeated_gi_full, significant_types_of_motifs,
# 								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
# 								  axis_text_colour=axis_text_colour)

false_positive_counts_threshold <- 0.02
p_value_threshold <- 0.05
p_value_column <- "enrichment_adj_p_values"

is_significant <- (observed_data[, "value"] > sum ( observed_data[, "value"] ) * false_positive_counts_threshold)  &
				  (p_values_data[, p_value_column]  < p_value_threshold)

significance_color <- sapply( is_significant, function(x) { return(boolean_to_colour(x)) } )


background_colour <- "white"   # "#7E4E99" # 
axis_text_colour <- "black"  # "#FFC000" # 

## Plot the plot

ggplot(random_data, aes(key, value)) + 
		geom_boxplot() + 
		geom_point( data=observed_data, aes(x=key, y=value  )   ## Add the observed values as additional points
					, color=significance_color, size=4 , alpha=0.5) 	+
		xlab( "Types of Motifs") 	+ 
		ylab( "Counts")  + 
		theme( plot.background = element_rect(fill = background_colour, color=background_colour), 
			   axis.text = element_text(colour = axis_text_colour),
			   axis.title = element_text(colour = axis_text_colour),
			   axis.text.x=element_text(face="italic"))
	
ggsave(file.path(figures_results_directory, "repeated_gi.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )


