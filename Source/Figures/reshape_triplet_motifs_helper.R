### Script: reshapte_triplet_motifs_helper.R
### Author: Ignatius Pang 
### Date: 13-7-2016
### Description: Helper functions for plotting results for poster presentation.
##                 * Compare observed and random values using box plot.
##				   * Unpivot the full results table from randomization analyses

# cd '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Poster'
# setwd ( '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Poster/')

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

print_box_plot_observed_vs_random <- function (observed_data, random_data, p_values_data, motifs_to_include = NULL, ordering=NULL, plot_type="boxplot", 
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
	
	## Fix missing motif types in observed data, according the ordering requested
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
	
	significance_color <- sapply ( is_significant, function(x) {  return(boolean_to_colour(x)) } )
	
	## Plot the plot
	
	if ( plot_type == "violin" ) {
	ggplot(random_data, aes(key, value)) + 
		geom_violin() + 
		geom_point( data=observed_data, aes(x=key, y=value  )   ## Add the observed values as additional points
					, color=significance_color, size=4 , alpha=0.5) 	+
		xlab( "Types of Motifs") 	+ 
		ylab( "Counts")  
	} else {
		
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
	}
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
	
	all_edge_types <- c('p', 'tu', 'td', 'ku', 'kd', 'na')
	
	all_edge_types_paper_style <- c('P', 'T1', 'T2', "K1", "K2", 'None')
	
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


