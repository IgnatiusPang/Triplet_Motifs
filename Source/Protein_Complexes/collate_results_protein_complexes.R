library(dplyr)
library(tidyr)
library(ggraptR)
library(ggplot2)
library(reshape2)

#####################################################################								 

### Graphic width and heigth for the boxplots of the randomization results. The observed value is shown in a dot. 
### This setting is for the poster presentation.

options <- commandArgs(trailingOnly = TRUE)
source( "./Common/parameters_file.R")

source_directory <- file.path( base_directory, "Source")
source_directory_common <- file.path( source_directory, "Common")
source( file.path(source_directory_common, "count_triplet_motifs_helper.R") )
source(  file.path(source_directory_common, "concatenate_result_files.R") )

results_directory                    <- file.path( results_directory, "Bootstrap_p_values/Protein_Complexes")
final_results_directory              <- file.path ( results_directory, "Final_Results")

create_directory_if_not_exists(final_results_directory)

################################################################################################################################################
# Counts of Triplet Motifs Costanzo et al. 2010 dataset

triplet_motifs_full_results_file     <- "full_results_table_benschop_protein_complex_in_tri_motifs.tab"
count_triplet_motifs_randomized_file <- "randomized_counts_benschop_protein_complex_in_tri_motifs.tab" 

triplet_motifs_file_pattern      	 <- "^randomized_counts_benschop_protein_complex_in_tri_motifs.tab"
observed_file_pattern        		 <-  "Job_20/observed_counts_benschop_protein_complex_in_tri_motifs.tab"
final_results_directory_pattern      <- "Final_Results"
number_of_randomized_samples         <- 2000
p_values 							 <- 0.05

collated_table <- collate_result_files_helper ( results_directory, triplet_motifs_full_results_file, count_triplet_motifs_randomized_file,
							  triplet_motifs_file_pattern, observed_file_pattern, final_results_directory_pattern,
							  number_of_randomized_samples, p_values)


