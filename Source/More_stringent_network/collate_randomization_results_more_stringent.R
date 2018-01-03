
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(purrr)

#####################################################################								 

## Load all the parameters 
options <- commandArgs(trailingOnly = TRUE)
source( "./Common/parameters_file.R")


### Graphic width and heigth for the boxplots of the randomization results. The observed value is shown in a dot. 
### This setting is for the poster presentation.
false_positive_rates <- 0.02

source_directory_common <- file.path( source_directory, "Common")
source( file.path(source_directory_common, "count_triplet_motifs_helper.R") )
source(  file.path(source_directory_common, "concatenate_result_files.R") )

################################################################################################################################################
# Counts of Triplet Motifs 3rd quartile

#folders_list <- Sys.glob ( "/home/ignatius/PostDoc/2016/Triplet_Motifs/Results/Bootstrap_p_values/More_stringent*")
# percentile <- c(100, seq ( 0.7, 0.95, 0.05)*100 )

percentile <- seq( 5, 100, 5) 

folders_list <- map_chr (percentile,  function(x) {
					 paste ( "/home/ignatius/PostDoc/2016/Triplet_Motifs/Results/Bootstrap_p_values/More_stringent_network_", x, sep="") })


my_collate_result <- function ( input_directory,  percentile  ) { 
	results_directory_percentile                    <- input_directory
	final_results_directory              <- file.path ( results_directory_percentile, "Final_Results")
	
	create_directory_if_not_exists(final_results_directory)
	
	triplet_motifs_file_pattern      	 <- "^randomized_results_table*"
	triplet_motifs_full_results_file     <- paste("full_results_triplet_motifs_more_stringent_", percentile ,"pc.tab", sep="")
	count_triplet_motifs_randomized_file <- paste("randomized_counts_triplet_motifs_more_stringent_", percentile ,"pc.tab",sep="" )
	
	job_number <- (percentile / 5*20)
	
	observed_file_pattern 	 	 <- paste( "Job_", job_number, "/results_count_triplet_motifs_costanzo.tab", sep="")
	final_results_directory_pattern  <- "Final_Results"
	number_of_randomized_samples         <- 2000
	p_values 						   	 <- 0.05
	
	final_results_table_count_triplet_motifs <- collate_result_files_helper ( results_directory_percentile, 
																			  triplet_motifs_full_results_file, count_triplet_motifs_randomized_file,
								  triplet_motifs_file_pattern, observed_file_pattern, final_results_directory_pattern, 
								  number_of_randomized_samples, p_values)

}

walk2( folders_list, percentile, my_collate_result)

################################################################################################################################################
