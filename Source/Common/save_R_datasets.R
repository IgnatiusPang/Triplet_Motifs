
library(RODBC)
library(igraph)
library(dplyr)
library(lazyeval)

# Author: Ignatius Pang
# Date: 17-5-2016
# Description: Transfer data from the SQL database into text file so it can be uploaded into R for data analysis without database connection.


# Rscript --vanilla save_R_datasets.R 

#################################################################################################################


options <- commandArgs(trailingOnly = TRUE)
source( "./Common/parameters_file.R")


source_directory_common  <- file.path( source_directory, "Common")

source( file.path(source_directory_common, "count_triplet_motifs_helper.R") )


#################################################################################################################

input_db     <- "sbi_triplet_motifs"  # The name of the SQL database for the data to be imported into
user_id      <- "ignatius"
user_passwd  <- "Sp33dBo@t"

conn_string  <- paste("SERVER=localhost;DRIVER={PostgreSQL};DATABASE=", input_db,
					 ";PORT=5432;UID=", user_id, ";PWD=", user_passwd ,";LowerCaseIdentifier=0", sep="")

channel      <- odbcDriverConnect( conn_string)

set_encoding <- sqlQuery(channel,"SET client_encoding='LATIN1'")

#################################################################################################################

gene_names_table <- sqlQuery ( channel, "select *
										   			 from  id_scer_gene_name 
										   			 where oln_id is not null 
										   				   and uniprot_acc is not null; ")

write.table (gene_names_table, file.path(data_directory, "gene_names_table.txt"), row.names = FALSE, sep="\t", quote=FALSE )


#################################################################################################################

hc_transcriptional_regulatory_network <- sqlQuery(channel, "select * from tf_high_confidence_transcriptional_regulation_network;", stringsAsFactors=FALSE)

write.table (hc_transcriptional_regulatory_network, file.path(data_directory, "tf_high_confidence_transcriptional_regulation_network.txt"),row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################

ppi_kinase_interaction_database_sharifpoor_2011 <- sqlQuery(channel, "select * from ppi_kinase_interaction_database_sharifpoor_2011;", stringsAsFactors=FALSE)

write.table (ppi_kinase_interaction_database_sharifpoor_2011, file.path(data_directory, "ppi_kinase_interaction_database_sharifpoor_2011.txt"),row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################

sbi_interactome <- sqlQuery(channel, "select * from results_ppi_sbi_hc_yeast_interactome_with_author_gene_name;", stringsAsFactors=FALSE)

write.table (sbi_interactome, file.path(data_directory, "results_ppi_sbi_hc_yeast_interactome_with_author_gene_name.txt"), row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################

costanzo_stringent <- sqlQuery(channel, "select * from genetic_interaction_costanzo_stringent_cutoff;", stringsAsFactors=FALSE)

write.table (costanzo_stringent, file.path(data_directory, "genetic_interaction_costanzo_stringent_cutoff.txt"), row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################

# Giaever et al. 2002 Functional profiling of the Saccharomyces cerevisiae genome. Nature. 2002 Jul 25;418(6896):387-91.
# pubmed ID: 12140549

# competitive fitness: decreased
# viable
# inviable

# 1100 essential genes (i.e. deletion of the gene leads to inviable phenotype)
essential_genes <- sqlQuery(channel, "select distinct oln_id from sgd_phenotype 
									  where reference ~ '12140549'
											and phenotype = 'inviable'  ; ", stringsAsFactors=FALSE)

write.table (essential_genes, file.path(data_directory, "essential_genes.txt"),row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################

# Sopko et al. 2006 (pubmed id: 16455487) List of proteins that are toxic if they are overexpressed. 
# 
# overexpressed_toxic_genes <- read.table ( "/media/babs/Systemsbiology/Yeast Systems Biology Database/Data/v1.02/Processed Data/Warehouse/Sopko 2006, 16455487/overexpressed_toxic_genes.csv", header=TRUE)
# 
# write.table (overexpressed_toxic_genes, file.path(data_directory, "overexpressed_toxic_genes.txt"), row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################
### Protein Complexes

benschop_protein_complexes <- sqlQuery(channel, "select distinct c_id as complex_id, oln_id from complex_ypd2complex_benschop;", stringsAsFactors=FALSE) 

write.table (benschop_protein_complexes, file.path(data_directory, "benschop_protein_complexes.txt"), row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################
### Periodically expressed genes

### Granovskaia 2010 et al. and Cyclebase
periodically_expressed_genes <- sqlQuery(channel, "select distinct oln_id
										from cell_cycle_peak_times_of_orfs; ") 

write.table (periodically_expressed_genes, file.path(data_directory, "periodically_expressed_genes.txt"), row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################

## Paralogs from OrthoMCL
orthomcl_paralogs <- sqlQuery(channel, "select distinct oln_id_a, oln_id_b from orthologue_group_orthomcl_scer_inparalogs;", stringsAsFactors=FALSE)
write.table (orthomcl_paralogs, file.path(data_directory, "orthomcl_paralogs.txt"), row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################

## ohnologs from SGD 

sgd_paralogs <- sqlQuery(channel, "select distinct oln_id_a, oln_id_b from sgd_paralogs;", stringsAsFactors=FALSE)
write.table (sgd_paralogs, file.path(data_directory, "sgd_paralogs.txt"), row.names = FALSE, sep="\t", quote=FALSE )

#################################################################################################################

tf_network 				<- hc_transcriptional_regulatory_network
kinase_network 			<- ppi_kinase_interaction_database_sharifpoor_2011
sbi_interactome 		<- sbi_interactome
gi_costanzo_stringent 	<- costanzo_stringent


#################################################################################################################

#### Actual observed number of triplet motifs

### Createa a table that include all kinase-substrate, protein-protein and transcription factor-gene interactions
### Do not contain genetic interactions

## All interactions in the kinase substrate network must be greater than this score
sharifpoor_kinase_substrate_network_threshold <- 2.5

## Thresholds for the Costanzo genetic interaction network
## Stringent negative genetic interactions must have score less than this threshold
costanzo_stringent_negative_gi_score_threshold    <- -0.12

## Stringent negative genetic interactions must have p-value less than this threshold
costanzo_stringent_negative_gi_p_value_threshold  <-  0.05	


costanzo_stringent_negative_gi_std_dev_threshold <- 0.08

##
tf_network_collated <- tbl_df(tf_network) %>%
						     clean_graph_table ( "regulator_oln_id", "target_oln_id", directed=TRUE) %>%
						collate_interactions_from_both_direction_tbl_df( "regulator_oln_id", "target_oln_id", 
																	   "transcription factor-target down", "transcription factor-target up",
																	   "td", "tu", 
																	   "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev")

write.table (tf_network_collated,  file.path(data_directory, "tf_network_collated.txt"), row.names = FALSE, sep="\t", quote=FALSE )


##
kinase_network_subset <- subset(kinase_network, score > sharifpoor_kinase_substrate_network_threshold)


write.table (kinase_network_subset,  file.path(data_directory, "kinase_network_subset.txt"), row.names = FALSE, sep="\t", quote=FALSE )

##
kinase_network_collated <- 	tbl_df(kinase_network_subset) %>%
     clean_graph_table ( "kinase_oln_id", "target_oln_id", directed=TRUE) %>%
	collate_interactions_from_both_direction_tbl_df( "kinase_oln_id", "target_oln_id", 
													 "kinase-substrate down", "kinase-substrate up",
													 "kd", "ku", 
													 "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev") 

write.table (kinase_network_collated,  file.path(data_directory, "kinase_network_collated.txt"), row.names = FALSE, sep="\t", quote=FALSE )


# Remove duplicated rows
sbi_interactome <-  dplyr::select ( sbi_interactome, one_of(c("oln_id_a", "oln_id_b")) ) %>% 
					dplyr::distinct()

sbi_interactome_collated <-  tbl_df(sbi_interactome) %>%
								clean_graph_table ( "oln_id_a", "oln_id_b", directed=FALSE) %>%
								collate_interactions_from_both_direction_tbl_df(  "oln_id_a", "oln_id_b", 
																			 "protein-protein", "protein-protein",
																			 "p", "p", 
																			 "oln_id_a", "oln_id_b", "interaction_type", "interaction_type_abbrev") 


write.table (sbi_interactome_collated,  file.path(data_directory, "sbi_interactome_collated.txt"), row.names = FALSE, sep="\t", quote=FALSE )


### All interactions except for genetic interactions
interactions_combined <- dplyr::bind_rows(tf_network_collated, kinase_network_collated, sbi_interactome_collated)

write.table (interactions_combined,  file.path(data_directory, "interactions_combined.txt"), row.names = FALSE, sep="\t", quote=FALSE )


### Applying the cut-off of abs ( std_dev ) < 0.08 does not change the results, but strictly speaking to leave them there.
filtered_costanzo_stringent_2010 <- tbl_df( subset(gi_costanzo_stringent, genetic_interaction_score < costanzo_stringent_negative_gi_score_threshold & 
											  	p_value < costanzo_stringent_negative_gi_p_value_threshold &
											  	abs ( std_dev ) < costanzo_stringent_negative_gi_std_dev_threshold
											  	))

write.table (filtered_costanzo_stringent_2010,  file.path(data_directory, "filtered_costanzo_stringent_2010.txt"), row.names = FALSE, sep="\t", quote=FALSE )


##
filtered_costanzo_stringent_2016 <- sqlQuery(channel, "select * from genetic_interaction_costanzo_stringent_cleaned_2016;", stringsAsFactors=FALSE)


write.table (filtered_costanzo_stringent_2016,  file.path(data_directory, "filtered_costanzo_stringent_2016.txt"), row.names = FALSE, sep="\t", quote=FALSE )


#################################################################################################################

save( gene_names_table, 
	  
	  tf_network, kinase_network, sbi_interactome, gi_costanzo_stringent, 
	  essential_genes,
	  # overexpressed_toxic_genes,
	  
	  benschop_protein_complexes,
	  periodically_expressed_genes,
	  orthomcl_paralogs, sgd_paralogs, 

	  tf_network_collated,
	  kinase_network_collated,
	  sbi_interactome_collated,
	  
	  kinase_network_subset,
	  interactions_combined,
	  filtered_costanzo_stringent_2010,
	  filtered_costanzo_stringent_2016,
	  sbi_interactome_collated,
	  
	  file=file.path(data_directory,  "network_data_library.Rdata") )

#################################################################################################################


