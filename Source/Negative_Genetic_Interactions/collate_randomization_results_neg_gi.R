
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

#####################################################################								 


## Load all the parameters 
options <- commandArgs(trailingOnly = TRUE)
source( "./Common/parameters_file.R")


### Graphic width and heigth for the boxplots of the randomization results. The observed value is shown in a dot. 
### This setting is for the poster presentation.
false_positive_rates <- 0.02

source_directory <- file.path( base_directory, "Source")
source_directory_common <- file.path( source_directory, "Common")
source( file.path(source_directory_common, "count_triplet_motifs_helper.R") )
source(  file.path(source_directory_common, "concatenate_result_files.R") )

results_directory                    <- file.path( base_directory, "Results/Bootstrap_p_values/Negative_Genetic_Interactions")
final_results_directory              <- file.path ( results_directory, "Final_Results")

################################################################################################################################################

# Counts of Triplet Motifs Costanzo et al. 2010 dataset
triplet_motifs_file_pattern      	 <- "^randomized_results_table*"
triplet_motifs_full_results_file     <- "full_results_triplet_motifs_costanzo_2016_collated.tab"
count_triplet_motifs_randomized_file <- "randomized_counts_triplet_motifs.tab" 
observed_file_pattern 	 	 <- "Job_21/results_count_triplet_motifs_costanzo.tab"
final_results_directory_pattern  <- "Final_Results"
number_of_randomized_samples         <- 2000
p_values 						   	 <- 0.05

final_results_table_count_triplet_motifs <- collate_result_files_helper ( results_directory, 
																		  triplet_motifs_full_results_file, count_triplet_motifs_randomized_file,
							  triplet_motifs_file_pattern, observed_file_pattern, final_results_directory_pattern, 
							  number_of_randomized_samples, p_values)

################################################################################################################################################

rows_to_include  <-  (final_results_table_count_triplet_motifs[, "observed_counts"] > sum ( final_results_table_count_triplet_motifs[, "observed_counts"] ) * false_positive_rates)  &
					(final_results_table_count_triplet_motifs[, "enrichment_adj_p_values"]  < p_values)

rows_to_include <- ifelse ( is.na(rows_to_include), FALSE, rows_to_include)

final_results_table_count_triplet_motifs[rows_to_include,]

################################################################################################################################################