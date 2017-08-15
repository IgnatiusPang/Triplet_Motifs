library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

#####################################################################								 

### Graphic width and heigth for the boxplots of the randomization results. The observed value is shown in a dot. 
### This setting is for the poster presentation.

base_directory <- "/home/ignatius/PostDoc/2016/Triplet_Motifs/"
figures_results_directory <- "/home/ignatius/PostDoc/2016/Triplet_Motifs/Results/Poster/Figures"
graphic_width <- 5
graphic_height <- 3

source_directory <- paste( base_directory, "Source/", sep="")
source_directory_common <- paste( source_directory, "Common/",  sep="")
source( paste(source_directory_common, "count_triplet_motifs_helper.R", sep="") )
source( paste(source_directory_common, "concatenate_result_files.R", sep="") )

results_directory       <- paste( base_directory, "Results/Bootstrap_p_values/Repeated_GI_in_Motifs_Freq_Dist/", sep="")
final_results_directory <- paste ( results_directory, "Final_Results/", sep="")

if (! file.exists(final_results_directory)){
		dir.create(final_results_directory)
}

################################################################################################################################################
# Counts of Triplet Motifs Costanzo et al. 2010 dataset

triplet_motifs_full_results_file     <- "full_results_counts_unique_gi_pairs_in_motifs.tab"
count_triplet_motifs_randomized_file <- "randomized_counts_unique_gi_pairs_in_motifs.tab" 

triplet_motifs_file_pattern      	 <- "^randomized_counts_unique_gi_pairs_in_motifs.tab"
observed_file_pattern        		 <-  "Job_20/observed_counts_unique_gi_pairs_in_motifs.tab"
final_results_directory_pattern      <- "Final_Results"
number_of_randomized_samples         <- 70197 
p_values 							 <- 0.05

randomized_table <- collate_randomization_result_files_into_one_table  (   results_directory, triplet_motifs_file_pattern,  final_results_directory_pattern, 
																  number_of_randomized_samples, id =TRUE   )

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

observed_data <- read.table ( file.path( results_directory, observed_file_pattern  ), header=TRUE)

observed_data <-  dplyr::filter(observed_data, motif_type!="others")

observed_data_cumsum <- observed_data %>%
					dplyr::group_by( motif_type) %>%
					dplyr::arrange( motif_type, counts) %>%
					dplyr::mutate( cumulative_count = cumsum(total_count))

randomized_table %>%
	dplyr::filter( motif_type != 'others') %>%
	ggplot (  aes ( counts, total_count, group=counts )) + 
	geom_boxplot() +
	#scale_y_log10() +
	geom_point( data= observed_data, aes( counts, total_count, color="red")) +
	facet_grid ( motif_type ~., scales="free" )

ggsave(file.path(figures_results_directory, "repeated_gi_num_motifs_per_gi.tiff"), plot=last_plot()) # , width=graphic_width, height=graphic_height 


randomized_table %>%
	dplyr::filter( motif_type != 'others') %>%
	dplyr::group_by(file_number, group, motif_type) %>%
	dplyr::mutate( cumulative_count = cumsum(total_count)) %>% 
	ggplot (  aes ( counts, cumulative_count, group=counts )) + 
	geom_boxplot() +
	#scale_y_log10() +
	geom_point( data= observed_data_cumsum, aes( counts, cumulative_count, color="red")) +
	facet_grid ( motif_type ~., scales="free" )	 
	
	
ggsave(file.path(figures_results_directory, "repeated_gi_num_motifs_per_gi_cumulative_sum.tiff"), plot=last_plot() ) # , width=graphic_width, height=graphic_height


