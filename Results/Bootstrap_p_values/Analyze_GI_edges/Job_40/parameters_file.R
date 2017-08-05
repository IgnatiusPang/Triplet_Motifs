### Script: count_triplet_motifs_negative_interactions.R
### Author: Ignatius Pang 
### Date: 20-5-2016
### Description: Count the observed number of triplet motifs. Randomize the respective networks and count the number triplet motifs, repeat randomization 1000 times.
### Obtain the bootstrap p-value

library(igraph)
library(dplyr)
library(lazyeval)
library(parallel)
library(tidyr)

sessionInfo()

#########################################################
# Source location

# psql sbi_triplet_motifs
# cd '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/'

# Command Line Parameters
#  script.R <compute location (katana|clive|local)>  <array job ID> <number of randomized trials>

#########################################################
### Parameters 

is_run_locally 			     <- TRUE
array_job_id 			     <- NULL 
random_number_seed 		     <- 1985
number_of_cores_to_use 		 <- 24
number_of_randomized_trials  <- 100
num_iteration_rewire_network <- NULL # Number of times each network is rewired before counting the triplet motifs
use_costanzo_2010_dataset 	 <- FALSE

data_directory 			<- "./"
source_directory 		<- "./"
results_directory 		<- "./"
source_directory_common <- "./"

if ( length(options) != 0 &
	! (  options[1] == 'katana' | options[1] == 'clive' | options[1] == 'locally'  )  ) {
		stop ( 'First command line parameter must either by katana, clive, or local') 
}

if ( length(options) != 0 
	 & ( options[1] == 'katana' | options[1] == 'clive')  ) {
	is_run_locally <- FALSE
	
	if ( options[1] == 'katana' ) {
		
		number_of_cores_to_use <- 16
		
	} else if (options[1] == 'clive' ) {
		
		number_of_cores_to_use <- 24
	}
}

print (paste ( "is_run_locally = ", is_run_locally) ) 


### Get array_job_id and set random number seed 
if ( is_run_locally == FALSE & length(options) >=2 ) {
	
	array_job_id <- as.numeric ( options[2]  )
	
	if ( ! is.null(array_job_id) ) {
		if ( is.numeric(array_job_id) ) { 
			random_number_seed <- random_number_seed + array_job_id 
		} else {
			stop ( "array_job_id is not numeric.")
		}
	}
}

## Deal with number of randomized trials
if ( is_run_locally == FALSE & length(options) >=3 ) {
	number_of_randomized_trials <- as.numeric ( options[3]  )
	
	if ( is.null(number_of_randomized_trials) | is.na(number_of_randomized_trials) ) {
		stop ( "number_of_randomized_trials is not numeric.")
	}
}

## Deal with parameters with the code are to be ran locally
if (is_run_locally == TRUE) {
	random_number_seed <- 1985
	number_of_cores_to_use <- 4
	number_of_randomized_trials <- 4
	num_iteration_rewire_network <- 0.001 # Number of times each network is rewired before counting the triplet motifs
	
	base_directory <- "/media/z3371724/PostDoc/2016/Triplet_Motifs/"

	data_directory <- paste( base_directory, "Data/Triplet_Motifs_R_data/", sep="")
	source_directory <- paste( base_directory, "Source/", sep="")
	source_directory_common <- paste( source_directory, "Common/",  sep="")
}

source( paste(source_directory_common, "count_triplet_motifs_helper.R", sep="") )

## Deal with random number generator: 
### Different documentation on how randomization will work when the processes are spared across different cores
# http://stackoverflow.com/questions/30456481/controlling-seeds-with-mclapply
# https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/RngStream.html
# http://www.informatik.uni-ulm.de/ni/staff/HKestler/Reisensburg2009/PDF/multicore.pdf
RNGkind("L'Ecuyer-CMRG")
set.seed(random_number_seed)

#########################################################

### Import datasets 

# Load the following tables containing data for biological networks:
# tf_network, 
# kinase_network, 
# sbi_interactome, 
# gi_costanzo_stringent, 
# gi_tf, 
# gi_kinase
load( paste(data_directory, "network_data_library.Rdata", sep=""))

if ( use_costanzo_2010_dataset == TRUE) {
	filtered_costanzo_stringent <- filtered_costanzo_stringent_2010
} else {
	filtered_costanzo_stringent <- filtered_costanzo_stringent_2016
}

#########################################################
