### Script: reshapte_triplet_motifs_results.R
### Author: Ignatius Pang 
### Date: 13-7-2016
### Description: Plotting results for poster presentation.
##                 * Compare observed value and frequency distribution of random values using box plots.
##				   * Compare observed value and the average values from randomization analyses using bar plots.  

# cd '/home/ignatius/PostDoc/2016/Triplet_Motifs/Source/Figures'

base_directory <- "/home/ignatius/PostDoc/2016/Triplet_Motifs"

setwd ( '/home/ignatius/PostDoc/2016/Triplet_Motifs/Source/Figures/')

library(dplyr)
library(tidyr)
# library(ggraptR)
library(ggplot2)
library(reshape2)

#####################################################################								 

source( file.path ( base_directory, 'Source/Figures/reshape_triplet_motifs_helper.R') )

computational_results_directory <- file.path ( base_directory, "Results/Bootstrap_p_values") 
figures_results_directory       <- file.path ( base_directory, "Results/Poster/Figures" ) 

plot_background_colour <- "white"   # "#7E4E99" # 
axis_text_colour <- "black"  # "#FFC000" # 


### Graphic width and heigth for the boxplots of the randomization results. The observed value is shown in a dot. 
### This setting is for the poster presentation.
graphic_width <- 5
graphic_height <- 3

#####################################################################								 

false_positive_rates <- 0.02

all_15_types_of_triplet_motifs <- c("pp" 
									,"tutu"
									,"pku"
									,"kuku"
									,"ptd"
									,"kdku"
									,"ptu"
									,"tdtd"
									,"pkd"
									,"tdtu"
									,"kdkd"
									,"tuku"
									,"tdku"
									,"tdkd"
									,"tukd")

################################################################################################
##### Using boxplot to show the random values and using a point to show the observed results  ##
################################################################################################

### Motifs that are significant. Sorted in decending order of the number of observed motifs.

negative_interactions_full <- read.table ( file= file.path( computational_results_directory, 
										   "Negative_Genetic_Interactions/Final_Results/full_results_triplet_motifs_costanzo_2016_collated.tab") )

rows_to_include  <-  (negative_interactions_full[, "observed_counts"] > sum ( negative_interactions_full[, "observed_counts"] ) * false_positive_rates)  &
					 (negative_interactions_full[, "enrichment_adj_p_values"]  < 0.05)

rows_to_include <- ifelse ( is.na(rows_to_include), FALSE, rows_to_include)

significant_types_of_motifs <- dplyr::arrange ( negative_interactions_full[rows_to_include, ], desc( observed_counts ) ) %>%
							   dplyr::select ( one_of ( c( "motif_type"))) %>% 
							   as.data.frame() %>%  
							   t() %>% 
							   as.vector()

## Convert motif names to paper format
significant_types_of_motifs <- convert_triplet_motifs_name_to_paper_style( significant_types_of_motifs ) 

# ################################################################################################

### Print the analysis of negative genetic interactions 
negative_interactions_full     <- read.table ( file= file.path( computational_results_directory, 
																"Negative_Genetic_Interactions/Final_Results/full_results_triplet_motifs_costanzo_2016_collated.tab") )

negative_interactions_full     <- dplyr::arrange ( negative_interactions_full, desc( observed_counts ) ) 
	
	
negative_interactions_random   <- read.table ( file=  file.path( computational_results_directory, 
																 "Negative_Genetic_Interactions/Final_Results/randomized_counts_triplet_motifs.tab"  ))

negative_interactions_observed <- read.table ( file= file.path( computational_results_directory, 
																"Negative_Genetic_Interactions/Job_21/results_count_triplet_motifs_costanzo.tab"), 
								               header=TRUE )

print_box_plot_observed_vs_random(negative_interactions_observed, negative_interactions_random, negative_interactions_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour)

ggsave(file.path(figures_results_directory, "significant_triplet_motifs_negative_gi.tiff"), plot=last_plot(), width=graphic_width, 
	   height=graphic_height ) 


# ################################################################################################

### Essential Genes, genetic interaction

essential_ab_full <- read.table ( file.path( computational_results_directory,
											 "Negative_GI_Essential/Final_Results/full_results_table_count_a_and_b.tab"), header =TRUE)

essential_ab_random   <- read.table ( file.path( computational_results_directory,
												 "Negative_GI_Essential/Final_Results/randomized_results_table_count_a_and_b.tab"), header =TRUE)

essential_ab_observed <- read.table ( file.path( computational_results_directory,
												 "Negative_GI_Essential/Job_20/observed_a_and_b_essential_gene_counts.tab"), header =TRUE)

print_box_plot_observed_vs_random(essential_ab_observed, essential_ab_random, essential_ab_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour)

ggsave(file.path(figures_results_directory, "essential_gene_ab.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )


# ### Essential Genes, 3rd protein

essential_c_full  <- read.table ( file.path( computational_results_directory,
												"Negative_GI_Essential/Final_Results/full_results_table_count_c.tab"), header =TRUE)

essential_c_random   <- read.table ( file.path( computational_results_directory,
											  "Negative_GI_Essential/Final_Results/randomized_results_table_count_c.tab"), header =TRUE)

essential_c_observed <- read.table ( file.path( computational_results_directory,
											  "Negative_GI_Essential/Job_20/observed_c_essential_gene_counts.tab"), header =TRUE)

print_box_plot_observed_vs_random(essential_c_observed, essential_c_random, essential_c_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour)

ggsave(file.path(figures_results_directory, "essential_gene_c.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )

################################################################################################
### Benschop Protein complexes
complexes_full   <- read.table ( file.path( computational_results_directory,
								 "Protein_Complexes/Final_Results/full_results_table_benschop_protein_complex_in_tri_motifs.tab"), header =TRUE)
complexes_random   <- read.table ( file.path( computational_results_directory,
								 "Protein_Complexes/Final_Results/randomized_counts_benschop_protein_complex_in_tri_motifs.tab"), header =TRUE)
complexes_observed <- read.table ( file.path( computational_results_directory,
								 "Protein_Complexes/Job_1/observed_counts_benschop_protein_complex_in_tri_motifs.tab"), header =TRUE)

print_box_plot_observed_vs_random(complexes_observed, complexes_random, complexes_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour)

ggsave(file.path(figures_results_directory, "benschop_protein_complexes.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )


################################################################################################
### Paralogs
orthomcl_paralogs_full   <- read.table ( file.path( computational_results_directory,
													  "Paralogs/Final_Results/full_results_table_orthomcl_paralogs_in_tri_motifs.tab"), header =TRUE)

orthomcl_paralogs_random   <- read.table ( file.path( computational_results_directory,
												   "Paralogs/Final_Results/randomized_counts_orthomcl_paralogs_in_tri_motifs.tab"), header =TRUE)
orthomcl_paralogs_observed <- read.table ( file.path( computational_results_directory,
												   "Paralogs/Job_1/observed_counts_orthomcl_paralogs_in_tri_motifs.tab"), header =TRUE)


print_box_plot_observed_vs_random(orthomcl_paralogs_observed, orthomcl_paralogs_random, orthomcl_paralogs_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour)

ggsave(file.path(figures_results_directory, "orthomcl_paralogs.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )


# 
# 
# orthomcl_paralogs_not_gi_full   <- read.table ( file.path( computational_results_directory, 
# 													"Paralogs/full_results_table_orthomcl_paralogs_not_gi_in_tri_motifs.tab"), header =TRUE)
# 
# orthomcl_paralogs_not_gi_random   <- read.table ( file.path( computational_results_directory, 
# 													  "Paralogs/randomized_counts_orthomcl_paralogs_not_gi_in_tri_motifs.tab"), header =TRUE)
# orthomcl_paralogs_not_gi_observed <- read.table ( file.path( computational_results_directory, 
# 													  "Paralogs/observed_counts_orthomcl_paralogs_not_gi_in_tri_motifs.tab"), header =TRUE)
# 
# 
# print_box_plot_observed_vs_random(orthomcl_paralogs_not_gi_observed, orthomcl_paralogs_not_gi_random, orthomcl_paralogs_not_gi_full,
# 								  significant_types_of_motifs, ordering=significant_types_of_motifs, plot_type="boxplot")
# 
# ggsave(file.path(figures_results_directory, "orthomcl_paralogs_not_in_gi.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height ) 
# 
# 
################################################################################################
### Paralogs not randomize GI
orthomcl_paralogs_full   <- read.table ( file.path( computational_results_directory,
													"Paralogs_Not_Randomize_GI/Final_Results/full_results_table_orthomcl_paralogs_in_tri_motifs.tab"), header =TRUE)

orthomcl_paralogs_random   <- read.table ( file.path( computational_results_directory,
													  "Paralogs_Not_Randomize_GI/Final_Results/randomized_counts_orthomcl_paralogs_in_tri_motifs.tab"), header =TRUE)
orthomcl_paralogs_observed <- read.table ( file.path( computational_results_directory,
													  "Paralogs_Not_Randomize_GI/Job_1/observed_counts_orthomcl_paralogs_in_tri_motifs.tab"), header =TRUE)


print_box_plot_observed_vs_random(orthomcl_paralogs_observed, orthomcl_paralogs_random, orthomcl_paralogs_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour, p_value_column="depletion_adj_p_values")

ggsave(file.path(figures_results_directory, "orthomcl_paralogs_not_randomize_gi.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )



################################################################################################
#### Unique GI

repeated_gi_full <- read.table ( file.path( computational_results_directory,
											"Unique_GI_in_Motifs/Final_Results/full_results_counts_unique_gi_pairs_in_motifs.tab"), header =TRUE)

repeated_gi_random   <- read.table ( file.path( computational_results_directory,
													  "Unique_GI_in_Motifs/Final_Results/randomized_counts_unique_gi_pairs_in_motifs.tab"), header =TRUE)
repeated_gi_observed <- read.table ( file.path( computational_results_directory,
													  "Unique_GI_in_Motifs/Job_1/observed_counts_unique_gi_pairs_in_motifs.tab"), header =TRUE)

print_box_plot_observed_vs_random(repeated_gi_observed, repeated_gi_random, repeated_gi_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour)

ggsave(file.path(figures_results_directory, "repeated_gi.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )

## Arguments used for testing the function
# p_values_data <- repeated_gi_full
# observed_data <- repeated_gi_observed
# random_data <- repeated_gi_random
# ordering <- NULL
# motifs_to_include <- significant_types_of_motifs


################################################################################################
#### Periodic gene expression

periodic_ab_full <- read.table ( file.path( computational_results_directory,
											"Cell_Cycle/Final_Results/full_results_periodic_gene_counts_a_and_b.tab"), header =TRUE)

periodic_ab_random   <- read.table ( file.path( computational_results_directory,
												"Cell_Cycle/Final_Results/randomized_results_periodic_gene_counts_a_and_b.tab"), header =TRUE)
periodic_ab_observed <- read.table ( file.path( computational_results_directory,
												"Cell_Cycle/Job_1/observed_a_and_b_periodic_gene_counts.tab"), header =TRUE)

print_box_plot_observed_vs_random(periodic_ab_observed, periodic_ab_random, periodic_ab_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour)

ggsave(file.path(figures_results_directory, "periodic_ab.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )


periodic_c_full <- read.table ( file.path( computational_results_directory,
											"Cell_Cycle/Final_Results/full_results_periodic_gene_counts_c.tab"), header =TRUE)

periodic_c_random   <- read.table ( file.path( computational_results_directory,
												"Cell_Cycle/Final_Results/randomized_results_periodic_gene_counts_c.tab"), header =TRUE)
periodic_c_observed <- read.table ( file.path( computational_results_directory,
												"Cell_Cycle/Job_1/observed_c_periodic_gene_counts.tab"), header =TRUE)

print_box_plot_observed_vs_random(periodic_c_observed, periodic_c_random, periodic_c_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour)

ggsave(file.path(figures_results_directory, "periodic_c.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )


################################################################################################
#### Gene Overexpression

overexpression_ab_full <- read.table ( file.path( computational_results_directory,
									  "Overexpressed_Toxic_Genes/Final_Results/full_results_table_overexpressed_toxic_gene_counts_a_and_b.tab"), header =TRUE)

overexpression_ab_random   <- read.table ( file.path( computational_results_directory,
										  "Overexpressed_Toxic_Genes/Final_Results/randomized_results_table_overexpressed_toxic_gene_counts_a_and_b.tab"), header =TRUE)

overexpression_ab_observed <- read.table ( file.path( computational_results_directory,
												"Overexpressed_Toxic_Genes/Job_1/observed_a_and_b_overexpressed_toxic_gene_counts.tab"), header =TRUE)

print_box_plot_observed_vs_random(overexpression_ab_observed, overexpression_ab_random, overexpression_ab_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour, 
								  axis_text_colour=axis_text_colour)

ggsave(file.path(figures_results_directory, "overexpression_ab.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )





overexpression_c_full <- read.table ( file.path( computational_results_directory,
										   "Overexpressed_Toxic_Genes/Final_Results/full_results_table_overexpressed_toxic_gene_counts_c.tab"), 
									       header =TRUE)

overexpression_c_random   <- read.table ( file.path( computational_results_directory,
										  "Overexpressed_Toxic_Genes/Final_Results/randomized_results_table_overexpressed_toxic_gene_counts_c.tab"),
										  header =TRUE)

overexpression_c_observed <- read.table ( file.path( computational_results_directory,
											   "Overexpressed_Toxic_Genes/Job_1/observed_c_overexpressed_toxic_gene_counts.tab"), header =TRUE)

print_box_plot_observed_vs_random(overexpression_c_observed, overexpression_c_random, overexpression_c_full, significant_types_of_motifs,
								  ordering=significant_types_of_motifs, plot_type="boxplot", background_colour=plot_background_colour )

ggsave(file.path(figures_results_directory, "overexpression_c.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height )


# ################################################################################################

## Randomly add or remove edges

# Need to run collate_random_edges_results.R first

### Directories
final_results_directory <- file.path( computational_results_directory, "Random_Edges/Final_Results") 

stat_test_result_random_edges <- read.table ( file=file.path( final_results_directory, "stat_test_result_random_edges.txt" ))

randomization_collated_counts <- read.table (  file=file.path ( final_results_directory, "randomization_collated_counts.txt"))


temp <- dplyr::select(stat_test_result_random_edges, one_of(c('parameter', 'motif_type', 'is_significant')))
temp$experiment <-1
temp$counts <-1
temp[, 'is_significant'] <- as.factor(temp[, 'is_significant'] )

randomization_collated_counts %>%
	dplyr::filter (  motif_type != 'total' ) %>%
	ggplot (  aes( experiment, counts) ) +
	geom_boxplot()  +	
	geom_rect(data = temp ,aes(  fill=is_significant ),
			  xmin = -Inf,xmax = Inf,
			  ymin = -Inf,ymax = Inf,alpha = 0.3) +
	scale_fill_manual( values=c('blue', 'red')) +
	facet_grid(   motif_type ~ parameter , scales="free") 


ggsave(file.path(figures_results_directory, "randomization_collated_results.tiff"), plot=last_plot(), width=graphic_width, height=graphic_height ) 

ggsave(file.path(figures_results_directory, "randomization_collated_results_15x15.tiff"), plot=last_plot(), width=15, height=15 ) 

# ################################################################################################


