#####################################################################								 

# Function: get_list_of_files
# In a directory containing many directories, with each directory containing results for one job, go down to the directory of each job and collate result files that matches a user supplided regexp pattern.

# Inputs:
# results_directory  - Directory that contain many sub-directories for job results
# file_pattern       - A string containing the regexp pattern

# Output:
# A list of paths. 
get_list_of_files <- function(results_directory, file_pattern  )  { 
	
	list_of_job_folders <-   list.files(path = results_directory)
	
	list_of_result_files <- c()
	
	for ( i in list_of_job_folders) {
		
		temp_file_list <- list.files(path = file.path( results_directory, i) )
		
		for ( j in temp_file_list ) {
			
			if (grepl(file_pattern, j) ) {
				
				one_file             <- file.path ( results_directory, i, "/", j )
				
				list_of_result_files <- c(list_of_result_files, one_file)  
			}
		}
	}
	
	return ( list_of_result_files)
}

#####################################################################								 

# Function: concatenate_table_from_file_list
# Given a list of files containing table with the same number of columns of the same types, concatenate the tables together into one.

# Input: 
# list_of_result_files - Complete file path for the result tables 
# id - add an id column indicating the file from which the file originates from 

# Output:
# result_table - A result table which has the same shape as the input tables.

concatenate_table_from_file_list <- function ( list_of_result_files, id=FALSE) {
	
	result_table <- c()
	file_counter <- 1
	
	for ( one_result_file in list_of_result_files ) {
		
		print( one_result_file)
		
		one_result_table  <- read.table( one_result_file, header=TRUE)
		
		if ( id == TRUE) {
			one_result_table <- dplyr::mutate(one_result_table, file_number= file_counter)
		}
		
		print ( length( one_result_table[,1]))
		
		if ( file_counter == 1) {
			result_table <- one_result_table
			
		} else {
			# Bind two tables together, filling in the columns that are missing. The two tables
			# will have the columns sorted in alphabetical / numerical order.
			result_table <-  bind_rows_of_two_tables  (result_table, one_result_table)
		}
		
		file_counter <- file_counter + 1
	}
	
	print ( length( result_table[,1]))
	
	
	return ( result_table)
}



#####################################################################

# Function: 
# Bind two tables together, filling in the columns that are missing. The two tables
# will have the columns sorted in alphabetical / numerical order.

# Input:
# result_table: The first table.
# one_result_table: The second table. The rows of the second table are to be concatenated to the end of the first one.

# Output:
# The two tables concatenated toghether. 
bind_rows_of_two_tables <- function (result_table, one_result_table) {
	
	
	## Add the missing columns if they are not present already.
	colnames_results_table <- colnames( result_table) 
	
	colnames_next_table <- colnames( one_result_table )
	
	add_to_next_table <- setdiff( colnames_results_table,  colnames_next_table)
	add_to_results_table <- setdiff( colnames_next_table, colnames_results_table)
	
	for ( i in add_to_results_table) {
		colnames_to_update <- c( colnames( result_table), i) 
		result_table <- cbind( result_table, rep( 0, length( result_table[,1])) )
		colnames( result_table) <- colnames_to_update
	}
	
	for ( j in add_to_next_table) {
		colnames_to_update <- c( colnames( one_result_table), j) 
		one_result_table <- cbind( one_result_table, rep( 0, length( one_result_table[,1])) )
		colnames( one_result_table) <- colnames_to_update
	}
	
	## Make sure the tables are sorted properly before binding them together
	result_table <- result_table[,sort( colnames( result_table) ) ]
	one_result_table <- one_result_table[, sort(colnames( one_result_table))]
	
	result_table <- rbind ( result_table, one_result_table)
	
	return ( result_table)
}



##################################################################################################################

# Function: collate_randomization_result_files_into_one_file
# 
# Inputs:
# triplet_motifs_full_results_file: name of the result file with the merged full bootstrapping results 
# count_triplet_motifs_randomized_file: name of the result file with the merged randomized counts
# triplet_motifs_file_pattern:  file pattern for the files with the randomized counts to merge
# observed_file_pattern:   file pattern with the observed counts to merge
# final_results_directory_pattern: file pattern for the directory where the final results are kept (i.e. do not merge files from this directory)
# number_of_randomized_samples: number of randomized samples to use 
# p_values : p-values threshold 

# Outpputs:
# A table with the results from many randomization runs collated. 

collate_randomization_result_files_into_one_table <- function (   raw_results_directory, triplet_motifs_file_pattern,  final_results_directory_pattern, 
										  number_of_randomized_samples = 2000, id = FALSE  ) { 
	
	
	triplet_motifs_list_of_result_files <- get_list_of_files (raw_results_directory, triplet_motifs_file_pattern  )  
	
	print( triplet_motifs_list_of_result_files)
	
	files_to_ignore <- grepl(final_results_directory_pattern, triplet_motifs_list_of_result_files)
	
	triplet_motifs_list_of_result_files <- triplet_motifs_list_of_result_files[ files_to_ignore == FALSE]
	
	files_to_use 						<- as.vector ( unlist ( sapply( triplet_motifs_list_of_result_files, function(x) { grepl( "Job", x) } )  ) )
	triplet_motifs_list_of_result_files <- triplet_motifs_list_of_result_files[files_to_use]
	count_triplet_motifs_randomized 	<- concatenate_table_from_file_list  ( triplet_motifs_list_of_result_files, id)
	
	if ( number_of_randomized_samples > length(count_triplet_motifs_randomized[,1] ) ) {
		
		stop_string <- paste( "number_of_randomized_samples (", number_of_randomized_samples, ") is greater than the number of rows available (",
							  length(count_triplet_motifs_randomized[,1] ), ").", sep="")
		
		stop ( stop_string)
	}
	
	count_triplet_motifs_randomized 	<- count_triplet_motifs_randomized[1:number_of_randomized_samples,]		

	return( count_triplet_motifs_randomized )	
}


# Function: collate_result_files_helper
# 
# Inputs:
# triplet_motifs_full_results_file: name of the result file with the merged full bootstrapping results 
# count_triplet_motifs_randomized_file: name of the result file with the merged randomized counts
# triplet_motifs_file_pattern:  file pattern for the files with the randomized counts to merge
# observed_file_pattern:   file pattern with the observed counts to merge
# final_results_directory_pattern: file pattern for the directory where the final results are kept (i.e. do not merge files from this directory)
# number_of_randomized_samples: number of randomized samples to use 
# p_values : p-values threshold 

# Outpputs:
# Table with the bootstrapping results 
## Table containing the following columns:
##    motif_type, observed_counts, 
##    enrichment_raw_p_values, enrichment_adj_p_values, 
##    depletion_raw_p_values, depletion_adj_p_values, 
##    mean, sd, median, number of randomized run, ratio

collate_result_files_helper <- function ( raw_results_directory, triplet_motifs_full_results_file, count_triplet_motifs_randomized_file,
										  triplet_motifs_file_pattern, observed_file_pattern, final_results_directory_pattern, 
										  number_of_randomized_samples = 2000 , p_values = 0.05) { 
	
	count_triplet_motifs_observed 	<- read.table ( file.path( raw_results_directory, observed_file_pattern ), 
													header=TRUE )
	
	count_triplet_motifs_randomized <- collate_randomization_result_files_into_one_table( raw_results_directory, 
																						  triplet_motifs_file_pattern,  
																						  final_results_directory_pattern, 
																						  number_of_randomized_samples = number_of_randomized_samples  )
	
	final_results_table_count_triplet_motifs <- get_full_results_table(count_triplet_motifs_observed, count_triplet_motifs_randomized, 
																	   p.value=p_values) 

	write.table ( count_triplet_motifs_randomized, file=file.path( raw_results_directory, final_results_directory_pattern, count_triplet_motifs_randomized_file ))
	write.table (final_results_table_count_triplet_motifs, file=file.path( raw_results_directory, final_results_directory_pattern, triplet_motifs_full_results_file ) )
	
	return( final_results_table_count_triplet_motifs)
}

################################################################################################################################################
