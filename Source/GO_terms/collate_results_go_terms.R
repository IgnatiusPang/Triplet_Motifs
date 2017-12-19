library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

#####################################################################								 

### Graphic width and heigth for the boxplots of the randomization results. The observed value is shown in a dot. 
### This setting is for the poster presentation.

options <- commandArgs(trailingOnly = TRUE)
source( "./Common/parameters_file.R")

source_directory_common <- file.path( source_directory, "Common")
source( file.path(source_directory_common, "count_triplet_motifs_helper.R") )
source(  file.path(source_directory_common, "concatenate_result_files.R") )

################################################################################################################################################
# Counts of Triplet Motifs, all three proteins share the same Gene Ontology Biological Process

results_directory_bp                    <- file.path( results_directory, "Bootstrap_p_values/GO_Terms_BP")
final_results_directory_bp              <- file.path ( results_directory_bp, "Final_Results")

create_directory_if_not_exists(final_results_directory_bp)

triplet_motifs_full_results_file     <- "full_results_table_go_bp_in_tri_motifs.tab"
count_triplet_motifs_randomized_file <- "randomized_counts_go_bp_in_tri_motifs.tab" 

triplet_motifs_file_pattern      	 <- "^randomized_counts_go_terms.*_in_tri_motifs.tab"
observed_file_pattern        		 <-  "Job_20/observed_counts_go_termsBP_in_tri_motifs.tab"
final_results_directory_pattern      <- "Final_Results"
number_of_randomized_samples         <- 2000
p_values 							 <- 0.05

collated_table <- collate_result_files_helper ( results_directory_bp, triplet_motifs_full_results_file, count_triplet_motifs_randomized_file,
							  triplet_motifs_file_pattern, observed_file_pattern, final_results_directory_pattern,
							  number_of_randomized_samples, p_values)



################################################################################################################################################
# Counts of Triplet Motifs, all three proteins share the same Gene Ontology Cellualar Compartment

results_directory_cc                    <- file.path( results_directory, "Bootstrap_p_values/GO_Terms_CC")
final_results_directory_cc              <- file.path ( results_directory_cc, "Final_Results")

create_directory_if_not_exists(final_results_directory_cc)

triplet_motifs_full_results_file     <- "full_results_table_go_cc_in_tri_motifs.tab"
count_triplet_motifs_randomized_file <- "randomized_counts_go_cc_in_tri_motifs.tab" 

triplet_motifs_file_pattern      	 <- "^randomized_counts_go_terms.*_in_tri_motifs.tab"
observed_file_pattern        		 <-  "Job_20/observed_counts_go_termsCC_in_tri_motifs.tab"
final_results_directory_pattern      <- "Final_Results"
number_of_randomized_samples         <- 2000
p_values 							 <- 0.05

collated_table <- collate_result_files_helper ( results_directory_cc, triplet_motifs_full_results_file, count_triplet_motifs_randomized_file,
												triplet_motifs_file_pattern, observed_file_pattern, final_results_directory_pattern,
												number_of_randomized_samples, p_values)

################################################################################################################################################
# Counts of Triplet Motifs, all three proteins share the same Gene Ontology Molecular Function

results_directory_mf                    <- file.path( results_directory, "Bootstrap_p_values/GO_Terms_MF")
final_results_directory_mf              <- file.path ( results_directory_mf, "Final_Results")

create_directory_if_not_exists(final_results_directory_mf)

triplet_motifs_full_results_file     <- "full_results_table_go_mf_in_tri_motifs.tab"
count_triplet_motifs_randomized_file <- "randomized_counts_go_mf_in_tri_motifs.tab" 

triplet_motifs_file_pattern      	 <- "^randomized_counts_go_terms.*_in_tri_motifs.tab"
observed_file_pattern        		 <-  "Job_20/observed_counts_go_termsMF_in_tri_motifs.tab"
final_results_directory_pattern      <- "Final_Results"
number_of_randomized_samples         <- 2000
p_values 							 <- 0.05

collated_table <- collate_result_files_helper ( results_directory_mf, triplet_motifs_full_results_file, count_triplet_motifs_randomized_file,
												triplet_motifs_file_pattern, observed_file_pattern, final_results_directory_pattern,
												number_of_randomized_samples, p_values)

################################################################################################################################################

