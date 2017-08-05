
### Script: test_random_edges.R
### Author: Ignatius Pang 
### Date: 20-5-2016
### Description: Test random additions or deletion of edges 

library(igraph)
library(dplyr)
library(lazyeval)
library(parallel)

sessionInfo()

#########################################################
# Source location

# psql sbi_triplet_motifs
# cd '/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/'
# 

#########################################################
### Parameters 
options <- commandArgs(trailingOnly = TRUE)

is_rm244 <- TRUE

if ( length(options) != 0 
     & options[1] == 'katana') {
  is_rm244 <- FALSE
}

print (is_rm244)

random_number_seed <- 1985
number_of_cores_to_use <- 16
number_of_randomized_trials <- 2000
num_iteration_rewire_network <- NULL # Number of times each network is rewired before counting the triplet motifs

data_directory <- "./"
source_directory <- "./"
results_directory <- "./"
source_directory_random_edges <- "./"

before_transpose_triplet_motif_counts <- "results_count_triplet_motifs_before_transpose_costanzo.tab"
output_observed_triplet_motif_counts <- "results_count_triplet_motifs_costanzo.tab" # "results_count_triplet_motifs_costanzo.tab"
randomized_results_table_file <- "randomized_results_table.tab"
output_enrichment_p_value <- "enrichment_p_value_fixed.tab"
output_depletion_p_value <- "depletion_p_value_fixed.tab"
output_full_results_table <- "full_results_table.tab"

if (is_rm244) {
  random_number_seed <- 1985
  number_of_cores_to_use <- 4
  number_of_randomized_trials <- 4
  num_iteration_rewire_network <- 0 # Number of times each network is rewired before counting the triplet motifs
  
  #base_directory <- "C:/Users/home/Desktop/Triplet_Motifs/"
   base_directory <- "/media/z3371724/PostDoc/2016/Triplet_Motifs/"
  
  data_directory <- paste( base_directory, "Data/Triplet_Motifs_R_data/", sep="")
  source_directory <- paste( base_directory, "Source/", sep="")
  results_directory <- paste( base_directory, "Results/Bootstrap_p_values/Random_Edges/", sep="")
  source_directory_random_edges <- paste( source_directory, "Random_Edges/", sep="")
}

# ## All interactions in the kinase substrate network must be greater than this score
# sharifpoor_kinase_substrate_network_threshold <- 2.5

# ## Thresholds for the Costanzo genetic interaction network
# ## Stringent negative genetic interactions must have score less than this threshold
# costanzo_stringent_negative_gi_score_threshold    <- -0.12

# ## Stringent negative genetic interactions must have p-value less than this threshold
# costanzo_stringent_negative_gi_p_value_threshold  <-  0.05	


source( paste(source_directory, "count_triplet_motifs_helper.R", sep="") )
source( paste (source_directory_random_edges, "random_edges_helper.R", sep="" ))

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


#########################################################

sbi_multiple_edges_no_rep <- get_multiple_edges  (sbi_interactome, "oln_id_a", "oln_id_b", duplicates=FALSE  )
	
sbi_multiple_edges_with_rep <- get_multiple_edges  (sbi_interactome, "oln_id_a", "oln_id_b", duplicates=TRUE  ) 
		
length ( sbi_multiple_edges_no_rep[,1])
length ( sbi_multiple_edges_with_rep[,1])



costanzo_multiple_edges_no_rep <- get_multiple_edges  (filtered_costanzo_stringent, "query_oln_id_edited", "array_oln_id_edited", 
													   directed=FALSE, duplicates=FALSE  )

costanzo_multiple_edges_with_rep <- get_multiple_edges  (filtered_costanzo_stringent, "query_oln_id_edited", "array_oln_id_edited", 
														 directed=FALSE, duplicates=TRUE  ) 

count ( filtered_costanzo_stringent)
count ( costanzo_multiple_edges_no_rep)
count ( costanzo_multiple_edges_with_rep)


#########################################################


sbi_interactome_edges_deleted <- randomly_remove_edges_from_network( sbi_interactome, 0.2)

length ( sbi_interactome_edges_deleted[,1])
length ( sbi_interactome[,1])*( 1-0.2)



kinase_network_subset_edges_deleted <- randomly_remove_edges_from_network( kinase_network_subset, 0.2)

length ( kinase_network_subset_edges_deleted[,1])
length ( kinase_network_subset[,1]) * ( 1 - 0.2)


#########################################################

#########################################################


sbi_interactome_edges_added<- randomly_add_edges_to_network (  sbi_interactome,  2, "oln_id_a", "oln_id_b", directed= FALSE,  
                                                               maximum_num_of_tries = 10 )

length ( sbi_interactome_edges_added[,1])
length ( sbi_interactome[,1]) * (2 + 1)



filtered_costanzo_stringent_cleaned <- clean_up_undirected_network( filtered_costanzo_stringent, "query_oln_id_edited", "array_oln_id_edited")

filtered_costanzo_stringent_edges_added <- randomly_add_edges_to_network (  filtered_costanzo_stringent_cleaned,  1, "query_oln_id_edited", "array_oln_id_edited", 
															   directed= FALSE,  
															   maximum_num_of_tries = 100 )

count ( filtered_costanzo_stringent_edges_added)
count ( filtered_costanzo_stringent_cleaned ) * (1 + 1)


#########################################################

kinase_network_subset_edges_added<- randomly_add_edges_to_network (  kinase_network_subset,  1, "kinase_oln_id", "target_oln_id", 
                                                                     directed= TRUE,  
                                                                     maximum_num_of_tries = 100 )

count (  kinase_network_subset_edges_added)[[1]]
count (  kinase_network_subset)[[1]]* (1 + 1)



tf_network_edges_added<- randomly_add_edges_to_network (  tf_network,  0.8, "regulator_oln_id", "target_oln_id",
																	 directed= TRUE,  
																	 maximum_num_of_tries = 100 )

count ( tf_network_edges_added)[[1]]
count ( tf_network)[[1]] * (1 + 0.8)


#########################################################

# List of all kinases
signalling_kinases <- unique( kinase_network[,"kinase_oln_id"])
signalling_targets <- unique( kinase_network[,"target_oln_id"])


# List of all transcription regulation factors
tf_regulators <- unique( tf_network[,"regulator_oln_id"])
tf_targets    <- unique( tf_network[,"target_oln_id"])


# Negative Genetic interactions
gi_query <- unique( as.vector ( t(  as.data.frame( filtered_costanzo_stringent[, "query_oln_id_edited"]  ) ) ) ) 
gi_array <- unique( as.vector ( t(  as.data.frame ( filtered_costanzo_stringent[, "array_oln_id_edited"] ) ) ) ) 

# > length ( unique ( gi_query,gi_array ))
# [1] 1645

# Protein-protein interactions
sbi_oln_id_a <- unique( sbi_interactome[, "oln_id_a"])
sbi_oln_id_b <- unique( sbi_interactome[, "oln_id_b"])

# All the genes in the network 
all_genes_in_network <- unique(c(signalling_kinases, signalling_targets, 
								 tf_regulators, tf_targets, 
								 gi_query, gi_array, 
								 sbi_oln_id_a, sbi_oln_id_b) ) 


#########################################################
### Testing at home
#  setdiff( triplet_motifs_costanzo_kinase_edges_added, triplet_motifs_costanzo) 
# 
#  setdiff( triplet_motifs_costanzo, triplet_motifs_costanzo_kinase_edges_added )
# 
# 
# tbl_df( setdiff ( kinase_network_subset_edges_added, kinase_network_subset) )
# tbl_df( setdiff ( kinase_network_subset, kinase_network_subset_edges_added) ) %>% arrange ( kinase_oln_id, target_oln_id)
# 
# as.data.frame( setdiff( interactions_combined_kinase_edges_added, interactions_combined )  )
# 
# as.data.frame( setdiff( interactions_combined, interactions_combined_kinase_edges_added )  )
# 

## Testing at work
# hello_my_friend <- filter ( triplet_motifs_costanzo, (type_ac== 'ku' & type_bc== 'kd') | (type_ac== 'kd' & type_bc== 'ku')  )
#hello_my_friend_2 <- filter ( triplet_motifs_costanzo, (type_ac== 'ku' & type_bc== 'kd') | (type_ac== 'kd' & type_bc== 'ku')  )
# filter ( hello_my_friend,   oln_id_b == 'YBR059C' & oln_id_a == 'YBL007C'& oln_id_c== 'YIL095W' )
# filter ( hello_my_friend_2, oln_id_a == 'YBR059C' & oln_id_b == 'YBL007C'& oln_id_c== 'YIL095W' )

#set_differences <- setdiff( hello_my_friend_2, hello_my_friend)

# filter ( interactions_combined, oln_id_a=="YBR059C" & oln_id_b == "YBL007C") 
# filter ( interactions_combined, oln_id_a=="YBL007C" & oln_id_b == "YIL095W") 
# 
# filter ( kinase_network_collated, oln_id_a=="YBR059C" & oln_id_b == "YBL007C") 
# filter ( kinase_network_collated, oln_id_a=="YBL007C" & oln_id_b == "YIL095W") 
# filter ( filtered_costanzo_stringent, (query_oln_id_edited=="YBR059C" & array_oln_id_edited == "YBL007C" )   | 
# 		 								(query_oln_id_edited== "YBL007C" & array_oln_id_edited ==  "YBR059C"  ) )
# 
# filter ( triplet_motifs_costanzo,   oln_id_a == 'YBR059C' & oln_id_b == 'YBL007C'& oln_id_c== 'YIL095W' )
# filter ( triplet_motifs_costanzo,   oln_id_b == 'YBR059C' & oln_id_a == 'YBL007C'& oln_id_c== 'YIL095W' )
