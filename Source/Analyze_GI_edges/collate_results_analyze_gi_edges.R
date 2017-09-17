library(dplyr)
library(tidyr)
#library(ggraptR)
library(ggplot2)
library(reshape2)

#####################################################################								 

### Graphic width and heigth for the boxplots of the randomization results. The observed value is shown in a dot. 
### This setting is for the poster presentation.

#  figure_file_suffix <- ".tiff"
 figure_file_suffix <- ".pdf"

## Set the parameters for this project
options <- commandArgs(trailingOnly = TRUE)
source( "./Common/parameters_file.R")

source_directory <- file.path( local_base_directory, "Source")
source_directory_common <- file.path( source_directory, "Common")
source_figures_directory <- file.path( source_directory, "Figures")
source( file.path(source_directory_common, "count_triplet_motifs_helper.R") )
source( file.path(source_directory_common, "concatenate_result_files.R") )
source( file.path(source_figures_directory, "paper_figures_helper.R") )

results_directory       <- file.path( local_base_directory, "Results/Bootstrap_p_values/Analyze_GI_edges")
final_results_directory <- file.path ( results_directory, "Final_Results")

create_directory_if_not_exists(final_results_directory)

################################################################################################################################################
### Helper functions 

fix_table_column_names <- function ( input_table ) { 
	
	type_ac <- "type_ac"
	type_bc <- "type_bc"
	motif_type <- "motif_type"
	mutate_dots <- ~paste(type_ac, type_bc, sep="")
	input_table <- mutate_(input_table, .dots=setNames(list(mutate_dots), c(motif_type) ) ) 
	
	colnames(input_table)[ colnames(input_table ) == "NA."  ] <- "NA"
	
	input_table <- dplyr::select ( input_table, 
								   one_of(c("motif_type", "kd", "ku", "p", "td", "tu",   "NA")))
	
	return( input_table)
}

################################################################################################################################################
# Counts of Triplet Motifs Costanzo et al. 2010 dataset

triplet_motifs_full_results_file     <- "full_results_table_analyze_gi_edge_in_tri_motifs.tab"
count_triplet_motifs_randomized_file <- "randomized_results_table.tab" 

triplet_motifs_file_pattern      	 <- "^randomized_results_table.tab"
observed_file_pattern        		 <- "Job_20/observed_results_table_file.tab"
final_results_directory_pattern      <- "Final_Results"
number_of_randomized_samples         <- 58673
p_values 							 <- 0.05
false_positive_rates 			     <- 0.02  # If the observed data is less than two percent of the total for that category, mark as insiginificant

number_to_use_for_bonferroni_adjustment <- 36

### Observed
count_triplet_motifs_observed 		<- read.table ( file.path( results_directory, observed_file_pattern), 
												header=TRUE )

count_triplet_motifs_observed <- fix_table_column_names(count_triplet_motifs_observed) 

count_triplet_motifs_observed_gathered <- tidyr::gather ( count_triplet_motifs_observed, "edge_type", "counts", 2:7) %>%
										  dplyr::mutate ( counts= ifelse (is.na(counts), 0 ,counts) ) 

### Randomized 
analyze_gi_edges_counts <- collate_randomization_result_files_into_one_table  (   results_directory, triplet_motifs_file_pattern,  
																				  final_results_directory_pattern, 
																		 			number_of_randomized_samples   ) 


analyze_gi_edges_counts <- fix_table_column_names(analyze_gi_edges_counts) 


## Adjustment for tuku 
adjustment_for_tuku <-  analyze_gi_edges_counts %>%
	dplyr::group_by ( motif_type) %>%
	dplyr::summarise( num_rows=n()) %>% 
	dplyr::filter( motif_type== "tuku") %>%
	dplyr::select( one_of('num_rows')) %>%
	as.data.frame() %>%
	as.vector() 

num_rows_to_add <- 4000 - adjustment_for_tuku[1]

rows_to_add <- data.frame( motif_type=rep( "tuku", num_rows_to_add), 
						   kd=rep( 0, num_rows_to_add), 
						   ku=rep( 0, num_rows_to_add),
						   p=rep( 0, num_rows_to_add), 
						   td=rep( 0, num_rows_to_add),
						   tu=rep( 0, num_rows_to_add), 
						   "NA"=rep( 0, num_rows_to_add))

colnames(rows_to_add)[ colnames(rows_to_add ) == "NA."  ] <- "NA"


analyze_gi_edges_counts <- rbind( analyze_gi_edges_counts, rows_to_add)


## tuku motif has much lower counts. 
# analyze_gi_edges_counts %>%
# 	group_by ( motif_type) %>%
# 	summarise( num_rows=n())


analyze_gi_edges_counts_gathered <- tidyr::gather ( analyze_gi_edges_counts, "edge_type", "counts", 2:7) %>%
									dplyr::mutate ( counts= ifelse (is.na(counts), 0 ,counts) ) 

### Compare the two to get a p-value 

joined_table <- dplyr::left_join ( analyze_gi_edges_counts_gathered, count_triplet_motifs_observed_gathered, 
							       by =c ( "motif_type" = "motif_type", "edge_type"= "edge_type")) %>%
				dplyr::rename( randomized_counts = counts.x, observed_counts = counts.y )  


### Replace NA with zeros

joined_table <- dplyr::mutate (joined_table, randomized_counts= ifelse (is.na(randomized_counts), 0 ,randomized_counts) ) %>%
				dplyr::mutate ( observed_counts= ifelse (is.na(observed_counts), 0 ,observed_counts) ) 
	

## Observed Count if greater than Randomized count, 1, else 0 
joined_table_summarized <- dplyr::mutate( joined_table, obs_greater_than_random = ifelse(observed_counts > randomized_counts, 0, 1 )) %>%
						   dplyr::group_by ( motif_type, edge_type, observed_counts ) %>%
						   dplyr::summarise( p_value=sum(obs_greater_than_random )/4000,   
												  mean=mean(randomized_counts), 
												  sd =sd( randomized_counts), 
												  median = median( randomized_counts), 
												  mad = mad(randomized_counts ) ) %>%
						   as.data.frame()

### Filter for the motif types that I need 

joined_table_summarized_filtered <- dplyr::filter (joined_table_summarized, motif_type %in% c( "pp", 'kuku', 'tutu', 'tdp', 'ptd','pku', 'pkd' )  ) 

analyze_gi_edges_counts_gathered <- dplyr::filter (analyze_gi_edges_counts_gathered, motif_type %in% 
												  	c( "pp", 'kuku', 'tutu', 'tdp', 'ptd', 'pku', 'pkd' )  ) 


#################################################################################################################################################
### Bonferroni adjustment
joined_table_summarized_filtered <- dplyr::mutate (joined_table_summarized_filtered,  p_value = ifelse( p_value *number_to_use_for_bonferroni_adjustment >1, 
																								 1, 
																								 p_value*number_to_use_for_bonferroni_adjustment) )

## Convert significance to colour 
sum_per_motif_type <- dplyr::group_by (joined_table_summarized_filtered, motif_type ) %>%
					  dplyr::summarise( total_per_motif_type=sum(observed_counts) )

joined_table_summarized_filtered <- dplyr::left_join ( joined_table_summarized_filtered, sum_per_motif_type, by="motif_type") %>%
									dplyr::mutate( is_significant =  ifelse( p_value < p_values & 
																		  	total_per_motif_type*false_positive_rates < observed_counts, 1, 0 )    )


joined_table_summarized_filtered <- dplyr::mutate ( joined_table_summarized_filtered, color=ifelse ( is_significant==1, 'red', 'blue'))

#################################################################################################################################################

#### Clean motif names
joined_table_summarized_filtered <- dplyr::mutate ( joined_table_summarized_filtered, edge_type=ifelse(edge_type=='NA', "na", edge_type))
analyze_gi_edges_counts_gathered <- dplyr::mutate ( analyze_gi_edges_counts_gathered, edge_type=ifelse(edge_type=='NA', "na", edge_type))

joined_table_summarized_filtered[,"motif_type"] <- convert_triplet_motifs_name_to_paper_style( 
																	joined_table_summarized_filtered[,"motif_type"] )

analyze_gi_edges_counts_gathered[,"motif_type"] <- convert_triplet_motifs_name_to_paper_style( 
																	analyze_gi_edges_counts_gathered[,"motif_type"] )



## Clean edge type names
joined_table_summarized_filtered[,"edge_type"] <- convert_edge_type_name_to_paper_style( 
													joined_table_summarized_filtered[,"edge_type"] )

analyze_gi_edges_counts_gathered[,"edge_type"] <- convert_edge_type_name_to_paper_style( 
													analyze_gi_edges_counts_gathered[,"edge_type"] )

### Order the data 
# motif_type_levels_order <- c('PP', 'TUTU', 'PKU', 'KUKU',  'TDP', 'PKD')

motif_type_levels_order <- c('PP', 'TUTU', 'TDP', 'PKU', 'KUKU', 'PKD')

edge_type_levels_order  <- c( 'None', 'P', 'K1', 'K2',  'T1', 'T2')

# motif_type_levels_order <- convert_triplet_motifs_name_to_paper_style( motif_type_levels_order)
# edge_type_levels_order  <- convert_edge_type_name_to_paper_style(edge_type_levels_order)

joined_table_summarized_filtered <- dplyr::mutate ( joined_table_summarized_filtered, edge_type=ifelse(edge_type=='NA', "na", edge_type))
analyze_gi_edges_counts_gathered <- dplyr::mutate ( analyze_gi_edges_counts_gathered, edge_type=ifelse(edge_type=='NA', "na", edge_type))

joined_table_summarized_filtered[, "motif_type"] <- factor ( joined_table_summarized_filtered[,"motif_type"], levels=motif_type_levels_order)
analyze_gi_edges_counts_gathered[, "motif_type"] <- factor ( analyze_gi_edges_counts_gathered[,"motif_type"], levels=motif_type_levels_order)

joined_table_summarized_filtered[, "edge_type"] <- factor ( joined_table_summarized_filtered[,"edge_type"], levels=edge_type_levels_order)
analyze_gi_edges_counts_gathered[, "edge_type"] <- factor ( analyze_gi_edges_counts_gathered[,"edge_type"], levels=edge_type_levels_order)

#################################################################################################################################################

write.table ( joined_table_summarized_filtered, file=file.path( final_results_directory, triplet_motifs_full_results_file), 
			  row.names=FALSE, quote=FALSE )

write.table ( analyze_gi_edges_counts_gathered, file=file.path( final_results_directory, count_triplet_motifs_randomized_file), 
			  row.names=FALSE, quote=FALSE )

#################################################################################################################################################


