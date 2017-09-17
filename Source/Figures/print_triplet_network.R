### Script: print_triplet_network.R
### Author: Ignatius Pang 
### Date: 13-7-2016
## Description: I want to visualize the network, which consists of all the triplet motifs. Only edges in triplet motifs are shown.

library(RCy3)
library(dplyr)

#######################################################################################################################

## Load the root base directory for the entire project
options <- commandArgs(trailingOnly = TRUE)
source( "./Common/parameters_file.R")

source( file.path ( source_directory_common, 'draw_triplet_motif_network_helper.R') ) 

#######################################################################################################################
# Read through the list of all triplet motifs 

triplet_motifs_costanzo <- readRDS( file=file.path( list_of_triplets_directory, "triplet_motifs_costanzo_2016.Rdata") ) 

triplet_motifs_costanzo <- as.data.frame ( triplet_motifs_costanzo )



gene_names <- read.table ( file.path(data_directory, "gene_names_table.txt"), sep="\t", header=TRUE, stringsAsFactors = FALSE, 
						   na.strings="", quote="")


## Results from Spatial Analysis of Functional Enrichment (SAFE tool)

enriched_functions <- read.table ( file.path( results_directory, "SAFE/safe_exported_reports_20161123.txt-node_properties_annotation-highest.txt"), sep="\t", header=TRUE, stringsAsFactors = FALSE, 
								   na.strings="", quote="")

#######################################################################################################################

cw <- draw_triplet_motif_network  ( triplet_motifs_costanzo, 'All triplets network' ) 

triplet_motif_network_display_style(cw) 

cw <- add_node_attributes(cw, gene_names , "oln_id")

cw <- add_node_attributes ( cw, enriched_functions, "Node.label.ORF")

empty_label_table <- gene_names %>% dplyr::mutate( empty_label = "") %>% dplyr::select ( one_of( c( 'oln_id', 'empty_label')))

cw <- add_node_attributes ( cw, empty_label_table, "oln_id")

setDefaultNodeSize( cw, new.size=0.1 )

list_of_nodes <- nodes(getGraph(cw))

# setNodeLabelDirect (cw, list_of_nodes, rep('', length(list_of_nodes)))

setNodeLabelRule(cw, 'empty_label')

displayGraph(cw)

current_example_dir <- file.path( results_directory, "SAFE")

create_directory_if_not_exists( current_example_dir) 

saveImage ( cw, file.name= file.path( current_example_dir, "Cloud_plot_edges_new"),
			image.type='pdf')

#######################################################################################################################
