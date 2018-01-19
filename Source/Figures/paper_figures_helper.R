### Script: paper_figures_helper.R
### Author: Ignatius Pang 
### Date: 13-7-2016
### Description: Helper functions for plotting results for the manuscript
##                 * Compare observed and random values using box plot.
##				   * Unpivot the full results table from randomization analyses



#####################################################################

## unpivot the columns motif_type, observed, and random
unpivot_tables <- function(input_table, treatment_group)  { 
	input_table_tidy <- input_table %>% 
		dplyr::rename( observed=observed_counts, random=mean)  %>%
		dplyr::select( one_of(c("motif_type", "observed", "random") )) %>%
		tidyr::gather(  "obs_or_rand", "value", -motif_type )
	
	input_table_tidy <- cbind( input_table_tidy, group=
							   	rep(treatment_group, length(input_table_tidy[,1])))
	
	return(input_table_tidy )
}

#####################################################################

## Compare observed and random values, significant p-value = red, insignificant p-value = blue

significance_to_colour <- function ( p_value, threshold) {
	
	
	if (is.null(p_value)   ) {
		
		return("blue")
	} 
	
	if (  is.na( p_value) | is.infinite(p_value) | is.nan(p_value) ) {
		
		return("blue")
		
	}

	if ( p_value >= threshold ) {
		
		return("blue")
		
	} else if ( p_value < threshold ) {
		
		return("red")
		
	} 
	
	return("blue")
	
}


## Compare observed and random values, significant p-value = red, insignificant p-value = blue
boolean_to_colour <- function ( x) {
	
	
	if (is.null(x)   ) {
		
		return("blue")
	} 
	
	if (  is.na( x) | is.infinite(x) | is.nan(x) ) {
		
		return("blue")
		
	}
	
	if ( x == TRUE ) {
		
		return("red")
		
	}
	
	
	
	if ( x == FALSE ) {
		
		return("blue")
		
	}

	
	return("blue")
	
}


# Function: print_box_plot_observed_vs_random_helper
# Description:
# 	Compare observed and random values using box plot.
# 	Significant observed data is shown as a red dot, Insignificant observed data is shown as a blue dot.
# 	Frequency distribution of many randomized datasets is shown as a box plot.
# 	Y-axis represents the frequency counts, 
#   X-axis represents eah type of triplet motif. 
# Input values:
#   observed_data: Data table of the observed values. Each column is one type of motif. Should only be one row representing observed value
#   random_data: Data table of the values from randomization runs. Each column is one type of motif. Each row represents the values from one randomized network.
#   p_values_data: Full results data table containing the observed and randomized results. Must contain column 'motif_type' and 'enrichment_adj_p_values' 
#   motifs_to_include: An array listing the types of triplet motifs to include in the plot.
#   ordering: If order is null, then sort by decreasing order of the observed frequency. Otherwise, rank the motif types from left to right by the ordered specified.
#   plot_type: 'boxplot' for Box-and-whisters plot, 'violin' for violin plot
#   background_colour: Backgound color of the plot
#   axis_text_colour: Colour of the text for the axis labels
#   false_positive_counts_threshold: if the observed count is lower than a proportion of the total observed counts, then it is considered a false positive
#   p_value_threshold: p-value below which the item is considered to be statistically significant
#   p_value_column: which column from the p_values_data column to obtain the p-value, default value = "enrichment_adj_p_values"
## Return:
# The tidied up observed_data
# The tidied up random_data
# A list called is_significant, showing which motif is significant
print_box_plot_observed_vs_random_helper <- function (observed_data, random_data, p_values_data, motifs_to_include = NULL, ordering=NULL, plot_type="boxplot", 
											   background_colour="white", axis_text_colour="black", false_positive_counts_threshold = 0.02, p_value_threshold=0.05,
											   p_value_column="enrichment_adj_p_values") {
	
	## Clean motif names
	motifs_to_include <- convert_triplet_motifs_name_to_paper_style( motifs_to_include ) 
	
	colnames( observed_data) <- convert_triplet_motifs_name_to_paper_style( colnames( observed_data) ) 
	
	colnames( random_data) <- convert_triplet_motifs_name_to_paper_style( colnames( random_data) ) 
	
	p_values_data[,"motif_type"] <- convert_triplet_motifs_name_to_paper_style( p_values_data[,"motif_type"] )
	
 	
	## Fix missing motif types, add what is missing in the observed data
	motif_types_to_be_added <- setdiff ( unique(colnames( random_data)   ),   
										 unique(colnames( observed_data) )  )
	
	if (length(motif_types_to_be_added) > 0 ) { 
		
		current_observed_data_colnames <- colnames( observed_data) 
		
		for ( i in 1:length(motif_types_to_be_added) ) { 
			observed_data <- cbind ( observed_data, c(0) )
		}
		
		colnames(observed_data ) <- c(current_observed_data_colnames, motif_types_to_be_added) 
	}
	
	## Fix missing motif types in observed data, according to the ordering requested
	if ( !is.null(ordering) & is.vector(ordering) ) { 
		motif_types_to_be_added <- setdiff ( ordering,   
											 unique(colnames( observed_data) )  )
	
	
			if (length(motif_types_to_be_added) > 0 ) { 
				
				current_observed_data_colnames <- colnames( observed_data) 
				
				for ( i in 1:length(motif_types_to_be_added) ) { 
					observed_data <- cbind ( observed_data, c(0) )
				}
				
				colnames(observed_data ) <- c(current_observed_data_colnames, motif_types_to_be_added) 
			}
	}

	## Fix missing motif types, add what is missing in the random data
	motif_types_to_be_added <- setdiff ( unique(colnames( observed_data) ),   
										   unique(colnames( random_data)   ) )
	
	if (length(motif_types_to_be_added) > 0 ) { 
			current_randomized_data_colnames <- colnames( random_data) 
			
			for ( i in 1:length(motif_types_to_be_added) ) { 
				random_data <- cbind ( random_data, rep(0, length( random_data[,1])) )
			}
			
			colnames(random_data ) <- c(current_randomized_data_colnames, motif_types_to_be_added) 
	}
	

	## Unpivot the data
	random_data <- gather( random_data)
	observed_data <- gather( observed_data)
	
	## Deal with the types of motifs to be included and how they are sorted in the boxplot
	if ( !is.null(motifs_to_include) 
		 & length(motifs_to_include) > 0 ) {
		
		random_data <- filter ( random_data, key %in% motifs_to_include)
		observed_data <- filter ( observed_data, key %in% motifs_to_include)
		p_values_data <- filter ( p_values_data, motif_type %in% motifs_to_include)
	}

	if ( is.null(ordering) ) {
		ordering  <- observed_data$key[order( observed_data$value, decreasing=TRUE)]
	} 
	
	random_data$key <- factor( random_data$key,  levels=ordering  ) 
	observed_data$key <- factor( observed_data$key, levels=ordering  ) 
		
	## re-arrange the data based on the ordering
	rownames( observed_data ) <- observed_data$key
	observed_data <- observed_data[ordering,]
	
	rownames( p_values_data) <- p_values_data[,"motif_type"]
	p_values_data <- p_values_data[ordering,]
		
	## Convert significance to colour 
	
	is_significant <- (observed_data[, "value"] > sum ( observed_data[, "value"] ) * false_positive_counts_threshold)  &
							(p_values_data[, p_value_column]  < p_value_threshold)
	
	names(is_significant) <-  p_values_data[,"motif_type"]
	
	return( list (  observed_data=observed_data,  
					random_data=random_data, 
					is_significant=is_significant  ))
}


# Function: print_box_plot_observed_vs_random
# Description:
# 	Compare observed and random values using box plot.
# 	Significant observed data is shown as a red dot, Insignificant observed data is shown as a blue dot.
# 	Frequency distribution of many randomized datasets is shown as a box plot.
# 	Y-axis represents the frequency counts, 
#   X-axis represents eah type of triplet motif. 
# Input values:
#   observed_data: Data table of the observed values. Each column is one type of motif. Should only be one row representing observed value
#   random_data: Data table of the values from randomization runs. Each column is one type of motif. Each row represents the values from one randomized network.
#   p_values_data: Full results data table containing the observed and randomized results. Must contain column 'motif_type' and 'enrichment_adj_p_values' 
#   motifs_to_include: An array listing the types of triplet motifs to include in the plot.
#   ordering: If order is null, then sort by decreasing order of the observed frequency. Otherwise, rank the motif types from left to right by the ordered specified.
#   plot_type: 'boxplot' for Box-and-whisters plot, 'violin' for violin plot
#   background_colour: Backgound color of the plot
#   axis_text_colour: Colour of the text for the axis labels
#   false_positive_counts_threshold: if the observed count is lower than a proportion of the total observed counts, then it is considered a false positive
#   p_value_threshold: p-value below which the item is considered to be statistically significant
#   p_value_column: which column from the p_values_data column to obtain the p-value, default value = "enrichment_adj_p_values"
#   log_y_axis: Use log y-axis if this equal to TRUE (default FALSE)

print_box_plot_observed_vs_random <- function (observed_data, random_data, p_values_data, motifs_to_include = NULL, ordering=NULL, plot_type="boxplot", 
											   background_colour="white", axis_text_colour="black", false_positive_counts_threshold = 0.02, 
											   p_value_threshold=0.05,
											   p_value_column="enrichment_adj_p_values", 
											   log_y_axis = FALSE) {

	

	
	# Add 1 to everything if using log y-axis
	if( log_y_axis==TRUE) {
		
		# Add 1 to observed data
		observed_data <- observed_data + 1
		
		# Change NA to zero
		random_data[is.na( random_data)] <- 0
		
		# Add 1 to randomized data
		random_data <- random_data + 1
	}
	
	tidy_data <- print_box_plot_observed_vs_random_helper(observed_data, random_data, p_values_data, motifs_to_include, 
														  ordering, plot_type, 
														  background_colour, axis_text_colour, 
														  false_positive_counts_threshold , p_value_threshold,
														  p_value_column) 
	
	is_significant <- tidy_data[["is_significant"]]
	observed_data <- tidy_data[["observed_data"]]
	random_data <- tidy_data[["random_data"]]
	
	
	significance_color <- sapply ( is_significant, function(x) {  return(boolean_to_colour(x)) } )
	

	
	## Plot the plot
	return_plot <- NULL 
	if ( plot_type == "violin" ) {
		return_plot <- ggplot(random_data, aes(key, value)) + 
		geom_violin() + 
		geom_point( data=observed_data, aes(x=key, y=value  )   ## Add the observed values as additional points
					, color=significance_color, size=4 , alpha=0.5) 	+
		xlab( "Types of Motifs") 	+ 
		ylab( "Counts")  
	} else {
		
		return_plot <- ggplot(random_data, aes(key, value)) + 
			geom_boxplot() + 
			geom_point( data=observed_data, aes(x=key, y=value  )   ## Add the observed values as additional points
						, color=significance_color, size=4 , alpha=0.5) 	+
			xlab( "Types of Motifs") 	+ 
			ylab( "Counts")  + 
			theme( plot.background = element_rect(fill = background_colour, color=background_colour), 
				   axis.text = element_text(colour = axis_text_colour),
				   axis.title = element_text(colour = axis_text_colour),
				   axis.text.x=element_text(face="italic"))
	}
	
	if( log_y_axis==TRUE) {
		return_plot <- return_plot + scale_y_log10()
	}
	
	
	return(return_plot )
}

#####################################################################

### Function: create_id_to_attribute_hash
### Description: Create a hash function that map keys to attributes. 

## Inputs: 
## keys: An array of key values
## attributes: An array of attribute values

## Output:
## An environment that act as a hash to convert keys to attributes. 

create_id_to_attribute_hash <- function(  keys, attributes) {
	
	keys <- as.character( as.vector(keys))
	attribute <- as.vector(attributes)
	
	hash <- new.env(hash = TRUE, parent = parent.frame())
	
	if ( length(keys) != length(attributes))  {
		warning('Length of keys != Length of attributes list.')
		return(1)
	}
	
	for ( i in 1:length(keys) )  {
		assign( keys[i], attributes[i], envir = hash)                       
	}
	
	return(hash)
}

##################################################################################################################

### Function: create_id_to_attribute_hash
### Description: Use a predefined hash dictionary to convert any Key to Attribute, return NA if key does not exists

## Inputs: 
## key: A key value
## hash: The hash dictionary that maps keys to attributes

## Output:
## A value that correspond to the query key value. 

convert_key_to_attribute <- function(key, hash) { 
	
	if ( exists(key, hash) ) {
		return ( get(key, hash)) 
	}
	else { 
		return (NA)
	}
}

##################################################################################################################

### Function: convert_keys_to_multiple_attributes
### Description: Use a predefined hash dictionary to convert multiple Keys to multiple attributes, 
###              return a specified value if key does not exists.

## Inputs: 
## keys: An array of key values
## hash: The hash dictionary that maps keys to attributes
## ifnotfound: A default value to return if a key does not exist

## Output:
## A array of attribute values that correspond to the key values. 


convert_keys_to_multiple_attributes <- function(keys, hash, ifnotfound=NA) { 
	
	return ( mget(keys, hash, ifnotfound=ifnotfound)) 
	
}


#####################################################################

### Function: convert_triplet_motifs_name_to_paper_style
# Inputs: 
# triplet_motif_names_array: Name of the triplet motifs in an array
# 
# Outputs:
# Convert the triplet motif names to the capitalized names used in the paper

# Testing Codes:
# test_table <- read.table ( '/media/z3371724/PostDoc/2016/Triplet_Motifs/Results/Bootstrap_p_values/Negative_Genetic_Interactions/Final_Results/full_results_triplet_motifs_costanzo_2016_collated.tab')
# 
# convert_triplet_motifs_name_to_paper_style ( test_table[,'motif_type'])


convert_triplet_motifs_name_to_paper_style <- function ( triplet_motif_names_array) {
	
	
	triplet_motif_names_array <- as.vector ( triplet_motif_names_array)
	
	all_15_types_of_triplet_motifs <- c("pp" 
										,"tutu"
										,"pku"
										,"kuku"
										,"tdp"
										,"kukd"
										,"tup"
										,"tdtd"
										,"pkd"
										,"tutd"
										,"kdkd"
										,"tuku"
										,"tdku"
										,"tdkd"
										,"tukd"
										### These are the ones that are already converted
										, "tdp"
										, "kdku"
										, "ptu"
										, "tdtu"
									
										, 'total'
										
										## Already capitalized
										, "PP"
										,"TUTU"
										,"PKU"
										,"KUKU"
										,"TDP"
										,"KDKU"
										,"PTU"
										,"TDTD"
										,"PKD"
										,"TDTU"
										,"KDKD"
										,"TUKU"
										,"TDKU"
										,"TDKD"
										,"TUKD"
										### These are the ones that are already converted
										, "TDP"
										, "KDKU"
										, "PTU"
										, "TDTU")
	
	all_15_types_of_triplet_motifs_paper_style <- c("PP" 
													,"TUTU"
													,"PKU"
													,"KUKU"
													,"TDP"
													,"KDKU"
													,"PTU"
													,"TDTD"
													,"PKD"
													,"TDTU"
													,"KDKD"
													,"TUKU"
													,"TDKU"
													,"TDKD"
													,"TUKD"
													### These are the ones that are already converted
													, "TDP"
													, "KDKU"
													, "PTU"
													, "TDTU"
													, 'total'
													
													## Already capitalized
													,"PP" 
													,"TUTU"
													,"PKU"
													,"KUKU"
													,"TDP"
													,"KDKU"
													,"PTU"
													,"TDTD"
													,"PKD"
													,"TDTU"
													,"KDKD"
													,"TUKU"
													,"TDKU"
													,"TDKD"
													,"TUKD"
													### These are the ones that are already converted
													, "TDP"
													, "KDKU"
													, "PTU"
													, "TDTU")
	
	hash <- create_id_to_attribute_hash(  all_15_types_of_triplet_motifs, 
										  all_15_types_of_triplet_motifs_paper_style) 
	
	
	updated_names <- convert_keys_to_multiple_attributes(triplet_motif_names_array, hash, ifnotfound=NA) 
	
	updated_names <- unlist (updated_names)
	
	
	if ( length( which ( is.na(updated_names ))) > 0 ) {
		
		stop_error_string <- 'convert_triplet_motifs_name_to_paper_style: List of input keys contain unmappable values: '
		
		unmappable_keys <- triplet_motif_names_array[is.na(updated_names )]
		
		stop_error_string <- paste( stop_error_string, unmappable_keys )
		
		stop ( stop_error_string)
	}
	
	return(updated_names )
}

#####################################################################

## Convert the name of edge types to the ones usable in the paper

convert_edge_type_name_to_paper_style <- function ( edge_type_names_array) {
	
	edge_type_names_array <- as.vector ( edge_type_names_array)
	
	all_edge_types <- c('p', 'tu', 'td', 'ku', 'kd', 'na', NA)
	
	all_edge_types_paper_style <- c('P', 'T1', 'T2', "K1", "K2", 'None', "None")
	
	hash <- create_id_to_attribute_hash(  all_edge_types, 
										  all_edge_types_paper_style) 
	
	updated_names <- convert_keys_to_multiple_attributes(edge_type_names_array, hash, ifnotfound=NA) 
	
	updated_names <- unlist (updated_names)
	
	if ( length( which ( is.na(updated_names ))) > 0 ) {
		
		stop_error_string <- 'convert_edge_type_name_to_paper_style: List of input keys contain unmappable values: '
		
		unmappable_keys <- edge_type_names_array[is.na(updated_names )]
		
		stop_error_string <- paste( stop_error_string, unmappable_keys )
		
		stop ( stop_error_string)
	}
	
	return(updated_names )
}




#####################################################################

# based on print_box_plot_observed_vs_random

clean_data_for_more_stringent <- function (observed_data, random_data, p_values_data, motifs_to_include = NULL, ordering=NULL, plot_type="boxplot", 
										   background_colour="white", axis_text_colour="black", false_positive_counts_threshold = 0.02, p_value_threshold=0.05,
										   p_value_column="enrichment_adj_p_values") {

		## Clean motif names
		motifs_to_include <- convert_triplet_motifs_name_to_paper_style( motifs_to_include ) 
		
		colnames( observed_data) <- convert_triplet_motifs_name_to_paper_style( colnames( observed_data) ) 
		
		colnames( random_data) <- convert_triplet_motifs_name_to_paper_style( colnames( random_data) ) 
		
		p_values_data[,"motif_type"] <- convert_triplet_motifs_name_to_paper_style( p_values_data[,"motif_type"] )
		
		
		## Fix missing motif types, add what is missing in the observed data
		motif_types_to_be_added <- setdiff ( unique(colnames( random_data)   ),   
											 unique(colnames( observed_data) )  )
		
		if (length(motif_types_to_be_added) > 0 ) { 
			
			current_observed_data_colnames <- colnames( observed_data) 
			
			for ( i in 1:length(motif_types_to_be_added) ) { 
				observed_data <- cbind ( observed_data, c(0) )
			}
			
			colnames(observed_data ) <- c(current_observed_data_colnames, motif_types_to_be_added) 
		}
		
		## Fix missing motif types in observed data, according to the ordering requested
		if ( !is.null(ordering) & is.vector(ordering) ) { 
			motif_types_to_be_added <- setdiff ( ordering,   
												 unique(colnames( observed_data) )  )
			
			
			if (length(motif_types_to_be_added) > 0 ) { 
				
				current_observed_data_colnames <- colnames( observed_data) 
				
				for ( i in 1:length(motif_types_to_be_added) ) { 
					observed_data <- cbind ( observed_data, c(0) )
				}
				
				colnames(observed_data ) <- c(current_observed_data_colnames, motif_types_to_be_added) 
			}
		}
		
		## Fix missing motif types, add what is missing in the random data
		motif_types_to_be_added <- setdiff ( unique(colnames( observed_data) ),   
											 unique(colnames( random_data)   ) )
		
		if (length(motif_types_to_be_added) > 0 ) { 
			current_randomized_data_colnames <- colnames( random_data) 
			
			for ( i in 1:length(motif_types_to_be_added) ) { 
				random_data <- cbind ( random_data, rep(0, length( random_data[,1])) )
			}
			
			colnames(random_data ) <- c(current_randomized_data_colnames, motif_types_to_be_added) 
		}
		
		
		## Unpivot the data
		random_data <- gather( random_data)
		observed_data <- gather( observed_data)

	
		return ( list( random_data = random_data, observed_data = observed_data))
}


############################



##########################################################



## Function: tidy_up_list_of_input_files

## Description: I use this to combine several plots together using faceting 
## Before I do this, I need to use this helper function to read the tables, and tidy up the significance value (e.g. apply 2% cutoff)
## Ideally use with map function to process lists of files. 
## Inputs:
# 
#  negative_interactions_full_file:     one pvalues_results_file
#  negative_interactions_random_file:	one	random_results_file
#  negative_interactions_observed_file:	one observed_results_file

# Outputs: 
#  A list containing the following elements: 
#   observed_results: a table with the following columns: motif_type, observed counts
#   random_results: a table with the following columns: motif_type, randomized counts
#   is_siginificant: a table with the following columns: motif_type, is_significant (TRUE or FALSE)

tidy_up_list_of_input_files <- function(negative_interactions_full_file, 
										negative_interactions_random_file,
										negative_interactions_observed_file
) { 
	### Print the analysis of negative genetic interactions 
	negative_interactions_full     <- read.table ( file= file.path( negative_interactions_full_file) )
	
	negative_interactions_full     <- dplyr::arrange ( negative_interactions_full, desc( observed_counts ) ) 
	
	
	negative_interactions_random   <- read.table ( file=  file.path( negative_interactions_random_file  ))
	
	negative_interactions_observed <- read.table ( file= file.path( negative_interactions_observed_file), 
												   header=TRUE )
	
	return_list_of_tidied_data <- print_box_plot_observed_vs_random_helper(negative_interactions_observed, negative_interactions_random, 
																		   negative_interactions_full, significant_types_of_motifs,
																		   ordering=significant_types_of_motifs, plot_type="boxplot", 
																		   background_colour=plot_background_colour, 
																		   axis_text_colour=axis_text_colour, 
																		   false_positive_counts_threshold = 0.02, p_value_threshold=0.05,
																		   p_value_column="enrichment_adj_p_values")
	
	return( return_list_of_tidied_data )
}


## Function: combine_graphs_using_faceting

## Description: I use this to combine several plots together using faceting 
## This saves me from using Inkscape to combine the several graphs into one multi-panel figure.

## Inputs:
# input_list: list of files with the following structure:
#    list(  pvalues_results_file_list, 
#			 random_results_file_list, 
#			 observed_results_file_list)

# list_of_facet_types:
#     List of categories for faceting with vertical stacking

# output_file_name: name of output file

# graphic_width: Width of the output graph

# graphic_height: Height of the output graph

# Outputs:
#  A list containing the following elements: 
#   observed_results: a table with the following columns: motif_type, observed counts, facet_type
#   random_results: a table with the following columns: motif_type, randomized counts, facet_type 
#   is_siginificant: a table with the following columns: motif_type, is_significant (TRUE or FALSE), facet_type    

combine_graphs_using_faceting <- function ( input_list, list_of_facet_types, sort_facet_decreasing =FALSE) {
	
	# Create the tidy up the list data tables
	list_of_tidied_data <- pmap( input_list, 
								 tidy_up_list_of_input_files)
	
	list_of_tidy_observed_data <- map(list_of_tidied_data, function(x) { x[["observed_data"]]} )
	list_of_tidy_random_data <- map(list_of_tidied_data, function(x) { x[["random_data"]]} )
	list_of_tidy_is_significant_vector <- map(list_of_tidied_data, function(x) { x[["is_significant"]]} )
	
	## work on observed_data
	list_of_tidy_observed_data_facet_type <- map2( list_of_tidy_observed_data, list_of_facet_types, 
												   function(x, y){  dplyr::mutate( x, facet_type=y )       }   )
	
	names( list_of_tidy_observed_data_facet_type) <- list_of_facet_types
	
	list_of_tidy_observed_data_facet_type <- list_of_tidy_observed_data_facet_type[ as.character(sort(list_of_facet_types, decreasing=sort_facet_decreasing)) ]
	
	## work on random_data
	list_of_tidy_random_data_facet_type  <- map2( list_of_tidy_random_data, list_of_facet_types, 
												  function(x, y){  dplyr::mutate( x, facet_type=y )       }   )
	
	names( list_of_tidy_random_data_facet_type ) <- list_of_facet_types
	
	list_of_tidy_random_data_facet_type <- list_of_tidy_random_data_facet_type[ as.character(sort(list_of_facet_types, decreasing=sort_facet_decreasing)) ]
	
	
	## work on is_significant
	
	transpose_is_significant_table_and_add_facet_types <- function(x, y){  
		
		transposed_data <- as.data.frame( x) %>% tibble::rownames_to_column()
		colnames( transposed_data) <- c( "motf_type", "is_significant") 
		
		transposed_data <- dplyr::mutate( transposed_data, facet_type=y )       
		
		return( transposed_data)
	} 
	
	list_of_tidy_is_significant_facet_type <- map2( list_of_tidy_is_significant_vector, list_of_facet_types, 
													transpose_is_significant_table_and_add_facet_types)
	
	names( list_of_tidy_is_significant_facet_type) <- list_of_facet_types
	
	list_of_tidy_is_significant_facet_type <- list_of_tidy_is_significant_facet_type[ as.character(sort(list_of_facet_types, decreasing=sort_facet_decreasing)) ]
	
	
	## Bind all rows from all tables together. Each table represents result from one percentile level
	observed_data_tidy    <- purrr::reduce (list_of_tidy_observed_data_facet_type, rbind)
	random_data_tidy      <- purrr::reduce (list_of_tidy_random_data_facet_type, rbind)
	is_significant_tidy   <- purrr::reduce (list_of_tidy_is_significant_facet_type, rbind)
	
	## Change facet as factor and relevel
	observed_data_tidy[,"facet_type"]  <- factor ( observed_data_tidy[,"facet_type"], levels = sort(list_of_facet_types, decreasing=sort_facet_decreasing ))
	random_data_tidy[,"facet_type"]    <- factor ( random_data_tidy[,"facet_type"], levels = sort(list_of_facet_types, decreasing=sort_facet_decreasing ))
	is_significant_tidy[,"facet_type"] <- factor ( is_significant_tidy[,"facet_type"], levels = sort(list_of_facet_types, decreasing=sort_facet_decreasing ))

	
	return( list ( observed_data=observed_data_tidy, 
				   random_data=random_data_tidy, 
				   is_significant=is_significant_tidy))	
}


## Input:
# observed_data, random_data, is_significant: These are outputsf from the 'combine_graphs_using_faceting' function
#   observed_results: a table with the following columns: motif_type, observed counts, facet_type
#   random_results: a table with the following columns: motif_type, randomized counts, facet_type 
#   is_siginificant: a table with the following columns: motif_type, is_significant (TRUE or FALSE), facet_type    
#   log_y_axis: Use log y-axis if this equal to TRUE (default FALSE)
# Outputs:
#    Graph saved as specified output file

motif_type_significance_box_plot_vertical_stack_faceting <- function( observed_data, random_data, is_significant, output_file_name, 
																	  graphic_width, graphic_height, log_y_axis=FALSE ) {
	
	significance_color <- sapply ( is_significant[,"is_significant"], function(x) {  return(boolean_to_colour(x)) } )
	
	background_colour <- "white"
	
	return_plot <- ggplot(random_data, aes(key, value)) + 
		geom_boxplot() + 
		geom_point( data=observed_data, aes(x=key, y=value  )   ## Add the observed values as additional points
					, color=significance_color, size=4 , alpha=0.5) 	+ # , size=4 , alpha=0.5
		xlab( "Types of Motifs") 	+ 
		ylab( "Counts")  + 
		facet_grid( facet_type ~., scales="free") +
		theme( plot.background = element_rect(fill = background_colour, color=background_colour), 
			   axis.text = element_text(colour = axis_text_colour),
			   axis.title = element_text(colour = axis_text_colour), 
			   axis.text.x=element_text(face="italic"))
	
	if( log_y_axis==TRUE) {
		return_plot <- return_plot + scale_y_log10()
	}
	
	ggsave(output_file_name, 
		   plot=return_plot, width=graphic_width, 
		   height=graphic_height ) 
}	


##########################################################



