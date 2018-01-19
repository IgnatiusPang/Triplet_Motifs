# Deciphering the network basis of negative genetic interactions in *Saccharomyces cerevisiae*with integrated biological networks and triplet motif analysis

# Quick summary 

* Analysis of Triplet Motifs consisting of genetic interactions, protein-protein interactions, transcription factor-gene target interactions and substrate-kinase interactions.

# Abstract

Negative genetic interactions in _Saccharomyces cerevisiae_ have been systematically screened to near-completeness, with >500,000 interactions identified. Nevertheless, the biological basis of these interactions remains poorly understood. To investigate this, we analyzed negative genetic interactions within an integrated biological network, being the union of proteinprotein, kinase-substrate, and transcription factor-target gene interactions. Network triplets, containing two genes / proteins that show negative genetic interaction and a third protein from the network, were then analyzed. Strikingly, just six out of 15 possible triplet motif types were present, as compared to randomized networks. These were in three clear groups: protein-protein interactions, signaling and regulatory triplets where the latter two showed no overlap. In the triplets, negative genetic interactions were associated with paralogs and ohnologs, however these were very rare. Negative genetic interactions among the six triplet motifs did however show strong dosage constraints, with genes being significantly associated with toxicity on overexpression and periodicity in the cell-cycle. Negative genetic interactions overlapped with other interaction types in 37% of cases; these were predominantly associated with protein complexes or signaling events. Finally, we highlight regions of ‘network vulnerability’ containing multiple negative genetic interactions; these could be targeted in fungal species for the regulation of cell growth. 

# Installations

## Software Required 
*R Statistical Computating Software:*  https://cran.r-project.org/
*R Studio, an Integrated Development Environment for R:* https://www.rstudio.com/ 
*Please install Java 8 first to use Cytoscape:* http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html
*Cytoscape - a software to visualize networks:* http://www.cytoscape.org/


## Please install these packages in R:

install.packages("ggpubr")
install.packages("igraph")
install.packages("knitr")  *  required for the 'kable' function for printing pretty table in html
* The parallel package is in the native R library, no installation required
install.packages("plyr") 
source("https://bioconductor.org/biocLite.R")
biocLite("RCy3") * The RCy3 package is reqiured for drawing network in Cytoscape using R scripts
install.packages("reshape2") 
install.packages("sqldf2") 
install.packages("svglite") 
install.packages("tidyverse") 


# Directories

* Negative_Genetic_Interactions
  + Analysis of triplet motifs containing one negative genetic interactions
  + Results Directory: Results/Bootstrap_p_values/Negative_Genetic_Interactions
  
* Random_Edges 
  + Random addition or removal of negative genetic interactions to test the robustness of the six overrepresented triplet motifs
  + Results Directory: Results/Bootstrap_p_values/Random_Edges
  
* More_stringent_network  
  + We repeated the triplet motif analysis but with data of increasingly stringent genetic interaction scores
  + Results Directory: Results/Bootstrap_p_values/More_stringent_network_X  (where X is a number within the range of 5 to 100)
  
* Paralogs
  + Pairs of proteins that participate in negative genetic interactions were examined to detect whether they were paralogs or ohnologs
  + Relevant Data Tables in 'Data/Triplet_Motifs_R_data/network_data_library.Rdata': orthomcl_paralogs, sgd_paralogs, 
  + Results Directory: Results/Bootstrap_p_values/Paralogs
  
* Analyze_GI_edges
  + Analysis of negative genetic interactions that overlap with some other types of interactions
  + Results Directory: Results/Bootstrap_p_values/Analyze_GI_edges

* Repeated_GI_in_Motifs
  + To find negative genetic interactions that are shared by two or more triplets, which were more frequent than by chance, their frequency was compared between observed and randomized networks
  + Results Directory: Results/Bootstrap_p_values/Repeated_GI_in_Motifs_Freq_Dist
  
* Negative_GI_Essential
  + To identify triplets that contain significant numbers of essential proteins
  + Relevant Data Table in 'Data/Triplet_Motifs_R_data/network_data_library.Rdata': essential_genes
  + Results Directory: Results/Bootstrap_p_values/Negative_GI_Essential
  
* Cell_Cycle
  + To test whether some triplet motifs are enriched for cell cycle-regulated genes 
  + Relevant Data Table in 'Data/Triplet_Motifs_R_data/network_data_library.Rdata': periodically_expressed_genes
  + Results Directory: Results/Bootstrap_p_values/Cell_Cycle
  
* Overexpressed_Toxic_Genes
  + To test whether some tiplet motifs were enriched for proteins that are toxi upon protein express
  + Relevant Data Table in 'Data/Triplet_Motifs_R_data/network_data_library.Rdata': overexpressed_toxic_genes
  + Results Directory: Results/Bootstrap_p_values/Overexpressed_Toxic_Genes
  
* GO_terms
  + To test whether the proteins that share negative genetic interaction are more likely to share the same GO term (GO slim terms)
  + Relevant Data Table in 'Data/Triplet_Motifs_R_data/network_data_library.Rdata': go_slim_mapping  
  + Results Directory: Results/Bootstrap_p_values/GO_Terms_XX (where XX is one of BP, CC, or MF. BP= GO Biological Process, CC = GO Cellular Compartment, MF = GO Molecular Function)
  
* Phenotype
  + To test whether the proteins that share negative genetic interaction are more likely to share the same phenotype
  + Relevant Data Table in 'Data/Triplet_Motifs_R_data/network_data_library.Rdata': phenotype_to_id_lookup_table
  + Results Directory: Results/Bootstrap_p_values/Phenotype  

* Common
  + Scripts and functions that are commonly shared between multiple scripts

* Figures
  + Scripts that list of all figures (and supplementary figuures) 
  + Results Directory: Results/Figures
  
* Supplementary_Files
  + Contain souce codes that were commonly used for the analysis
  + Results Directory: Results/Supplementary_Files
  
* Examples 
  + Some examples for analysis of triplet motifs

# Stucture of Analysis and Results Directory
Most of the data analysis directories contain the following script

* Calculation (R script) - The R script to perform the calculation. It uses the mclapply function to distribute randomisation jobs to multiple cores. 
* Run script (bash PBS script) - The Bash script to run the R script on the UNSW Katana cluster. Information on the Katana cluster (https://www.hpc.science.unsw.edu.au/cluster/katana). The cluster is managed by the PBS batch script system. 
* Collate script - Once I've got all the results from running different jobs on different cores and compute notes, I copy all the results to the results directory. There will be results from a number of job. This script collates all the results into one file. 
* The results directory often contains the results from multiple PBS jobs. Each job has its own directory (e.g. Job_X, where X is an integer). 
* The collated results are often saved in the directory named 'Final_Results'. This directory contains the full results table containing the observed counts for each triplet motif and the expected count from randomized networks. This directory also contains a file collating all the results from many randomisation tests.


# Version 
* Version 1.0 
* Project start 20th April 2016


* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)