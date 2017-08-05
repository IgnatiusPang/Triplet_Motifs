library(tidyr)
library(dplyr)
library(ggplot2)
library(lazyeval)


#####################################################################################################################################

collate_different_datasets <- function( results_directory, list_of_file_names, strsplit_pattern = '_', position_of_parameter = 6, expt_label) { 
	collated_table <- NULL
	
	prev_parameter <- -1000 # Dummy value 
	expt_number_counter <- 0 
	
	for ( i in 1:length ( list_of_file_names)) { 
		
		current_file_name <- list_of_file_names[i]
		
		file_name_splitted <- strsplit (current_file_name, strsplit_pattern )
		
		parameter_to_capture <- as.numeric( file_name_splitted[[1]][position_of_parameter] )
		
		if ( prev_parameter != parameter_to_capture ) {
			
			prev_parameter <- parameter_to_capture
			expt_number_counter <- 0
		}
		
		one_data_table <- read.table ( paste( results_directory, current_file_name, sep=""), header=TRUE)
		
		motif_types <- colnames( one_data_table)
		
		row_sums <- as.data.frame( rowSums( one_data_table, na.rm=TRUE) )
		colnames(row_sums) <- 'total'
		
		## Continue with the experiment number if the network size parameter is the same, otherwise reset to zero
		expt_number <- as.data.frame( (1+expt_number_counter):(length(one_data_table[,1] )+expt_number_counter) )
		
		expt_number_counter <- expt_number_counter + length ( one_data_table[,1])
	
		colnames( expt_number) <- 'expt_number'
		
		# print ( "hello")
		# print( dim(row_sums))
		# print (dim( expt_number))		
		# print ( dim(one_data_table ))
		
		one_data_table <- cbind( one_data_table, total=row_sums, number=expt_number)
		
		columns_to_get <- which ( colnames(one_data_table ) %in%  c( motif_types, 'total'))
		
		one_data_table_gathered <- tidyr::gather ( one_data_table, 'motif_type', 'counts',
												   columns_to_get)
		
		
		parameter_column_data <- rep( parameter_to_capture, length ( one_data_table_gathered[,1]) ) 
		expt_label_column_data <- rep( expt_label, length ( one_data_table_gathered[,1]) ) 
		
		one_data_table_gathered <- cbind( experiment= expt_label_column_data,
										  parameter = parameter_column_data,
										  one_data_table_gathered ) 
		
		if ( i == 1) {
			
			collated_table <- 	one_data_table_gathered
		} else {
			
			collated_table <- rbind ( collated_table, one_data_table_gathered)
			
		}
		
		
	}
	
	
	return (collated_table) 
}

#####################################################################################################################################



## Perform statistical testing

random_edges_statstical_testing <- function( randomization_collated_counts, 
											 counts_percentage_threshold = 0.02, p_value_threshold= 0.05, 
											 num_of_motifs=15, num_random_trials=2000, 
											 experiment_label_observed='Obs.',
											 experiment_label_random='Rnd.') { 
	
	parameters <- unique( randomization_collated_counts[, "parameter"] )
	
	motif_types <- setdiff ( unique( randomization_collated_counts[, "motif_type"] ), 'total' )
	
	result_parameters  <- c()
	result_motif_types <- c()
	result_p_values    <- c()
	results_95_quantile   <- c()
	results_prop_significant <- c()
	results_is_significant <- c()
	results_mean_before <- c()
	results_mean_after <- c()
	results_sd_before <- c()
	results_sd_after <- c()	
	
	for ( param in parameters ) {
		
		filter_condition_total_before <- interp( quote ( parameter==x &  experiment==y & motif_type=='total'), 
												 x=param, y= experiment_label_observed ) 
		
		total_before <- dplyr::filter( randomization_collated_counts, filter_condition_total_before  )	%>%
			dplyr::select(one_of ( c('expt_number', 'counts')))  %>%
			dplyr::rename( total_counts = counts)
		
		# print ( paste ( "dim(total_before) =", dim(total_before)))
		# print ( length(unique( total_before[,'expt_number'])) )
		
		if ( length(unique( total_before[,'expt_number'])) != length( total_before[,'expt_number'])) {
			stop( 'random_edges_statstical_testing: problem with unique key expt_number')
		}
		
		for ( motif_type in motif_types ) {
			
			
			filter_condition_before <- interp( quote ( parameter==x & motif_type==y & experiment==z), 
											   x=param, y=motif_type, z=experiment_label_observed) 
			
			before_table <- dplyr::filter( randomization_collated_counts, filter_condition_before  )
			
			# print ( paste( "dim(before_table) =", dim(before_table)) ) 
			
			## Counts must be at least 2% of total, otherwise reject
			####
			before_table_to_merge <- before_table %>%
				dplyr::select(one_of ( c('expt_number', 'counts')))  %>%
				dplyr::rename ( motif_counts = counts)
			
			# print ( "Before inner join")
			# print ( paste( "dim(before_table_to_merge) =", dim(before_table_to_merge)) ) 
			# print ( head( before_table_to_merge))
			
			before_table_to_merge <- inner_join ( before_table_to_merge, total_before, by=c( 'expt_number' = 'expt_number'))
			
			# print ( "After inner join")
			# print ( paste( "dim(before_table_to_merge) =", dim(before_table_to_merge)) ) 
			# print (head( before_table_to_merge))
			
			is_enough_frequency <- before_table_to_merge[,"motif_counts"] > before_table_to_merge[, "total_counts"]*counts_percentage_threshold
			####
			
			before <- before_table %>%
				dplyr::select(one_of ( c('counts'))) %>% 
				t() %>%
				as.vector()
			
			
			if ( length( before ) > 0 ){
				
				before[is.na(before)] <- 0 
			} else {
				
				before <- rep( 0, num_random_trials)
			}
			
			# print ( head(before))
			
			
			filter_condition_after <- interp( quote ( parameter==x & motif_type==y & experiment==z), 
											  x=param, y=motif_type, z=experiment_label_random) 
			
			after  <- dplyr::filter( randomization_collated_counts, filter_condition_after )	%>%
				dplyr::select(one_of ( c('counts'))) %>% 
				t() %>%
				as.vector()
			
			if ( length( after ) > 0 ){
				
				after[is.na(after)] <- 0 
			} else {
				
				after <- rep( 0, num_random_trials)
			}
			
			# print ( head(after))
			
			wilcox_test_result <-	wilcox.test ( before, after, 
												alternative="greater")$p.value * num_of_motifs # There are 15 types of triplet motifs 
			
			quantile_95_perc <- quantile(after, 1-(p_value_threshold/ num_of_motifs), na.rm=TRUE )[1]
			
			#### Check is above p-value threshold
			which_above_random_quartile <- which ( before > quantile_95_perc  )
			
			prop_above_random_quantile <-  (1 - length ( which_above_random_quartile  )/length( before) ) 
			
			#### Check is above p-value threshold and has enough observed counts
			accepted  <- which ( (before > quantile_95_perc) & is_enough_frequency )
			
			# print ( paste ( "before =", before, "quantile_95_perc =", quantile_95_perc, "is_enough_frequency =", is_enough_frequency))
			
			prop_siginifcant <- (1 - length ( accepted  )/length( before) ) 
			
			# print ( paste ( "length(is_enough_frequency)=", length(is_enough_frequency), 
			# 				"length(before) =", length(before), "length(accepted) =", length(accepted), 
			# 				"calc =", (1 - length ( accepted  )/length( before) ) )	)
			
			test_significance <- ifelse ( prop_siginifcant < p_value_threshold , 1, 0 )
			
			#### 
			results_mean_before <- c( results_mean_before, mean( before,na.rm = FALSE) ) 
			results_mean_after  <- c( results_mean_after, mean(after, na.rm = FALSE) )
			results_sd_before   <- c( results_sd_before, sd( before,na.rm = FALSE) )
			results_sd_after    <- c( results_sd_after, sd(after, na.rm = FALSE) )	
			
			results_95_quantile        <- c( results_95_quantile, prop_above_random_quantile )
			results_prop_significant   <- c( results_prop_significant, prop_siginifcant )
			results_is_significant     <- c( results_is_significant, test_significance )
			
			result_p_values     <- c( result_p_values, wilcox_test_result)
			result_parameters   <- c( result_parameters, param)
			result_motif_types  <- c( result_motif_types, motif_type)
		}
	}
	
	results_table <- data.frame ( parameter=result_parameters, 
								  motif_type=result_motif_types,
								  wilcox_adj_p_value=result_p_values, 
								  boot_adj_p_value=results_95_quantile,
								  boot_n_counts_adj_p_value= results_prop_significant,
								  is_significant= results_is_significant,
								  mean_before= results_mean_before, 
								  sd_before= results_sd_before,
								  mean_after= results_mean_after,
								  sd_after= results_sd_after)
	
	return( results_table)
	
}

#####################################################################################################################################