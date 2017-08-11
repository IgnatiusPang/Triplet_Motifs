# Count the number of each type of triplet motifs that each negative genetic interaction are involved in 

#########################################################
### Parameters 
options <- commandArgs(trailingOnly = TRUE)


### Source the parameters file here 
if ( length(options) != 0  )  { 
	if ( options[1] == 'katana' | options[1] == 'clive')   {
		source( "./parameters_file.R" )
		
	} else {
		source( "/media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Common/parameters_file.R")
	}
}

### Local Parameters
if (is_run_locally) {
	
	results_directory <- "/media/z3371724/PostDoc/2016/Triplet_Motifs/Results/Bootstrap_p_values/Unique_GI_in_Motifs/"
}

#########################################################

# triplet_motifs_costanzo

mutate_dots <- ~ifelse( type_ac >= type_bc, paste(type_ac, type_bc, sep=""), paste(type_bc, type_ac, sep="") ) 

triplet_motifs_costanzo_edited <- mutate_(triplet_motifs_costanzo, .dots=setNames(list(mutate_dots), c("motif_type") ) )

## Check all triplet motifs have oln_id_a >= oln_id_b
count(triplet_motifs_costanzo_edited) == filter ( triplet_motifs_costanzo_edited, oln_id_a >= oln_id_b) %>% count()

count_motif_per_negative_gi <- triplet_motifs_costanzo_edited %>%
									group_by ( oln_id_a, oln_id_b, motif_type) %>%
									summarise( count=n()) 

count_motif_per_negative_gi_spread <- count_motif_per_negative_gi %>% 
										spread( motif_type, count)
	

#########################################################






