# Deciphering the network basis of negative genetic interactions in *Saccharomyces cerevisiae*with integrated biological networks and triplet motif analysis

# Quick summary 

* Analysis of Triplet Motifs consisting of genetic interactions, protein-protein interactions, transcription factor-gene target interactions and substrate-kinase interactions.

# Abstract

Negative genetic interactions in _Saccharomyces cerevisiae_ have been systematically screened to near-completeness, with >500,000 interactions identified. Nevertheless, the biological basis of these interactions remains poorly understood. To investigate this, we analyzed negative genetic interactions within an integrated biological network, being the union of proteinprotein, kinase-substrate, and transcription factor-target gene interactions. Network triplets, containing two genes / proteins that show negative genetic interaction and a third protein from the network, were then analyzed. Strikingly, just six out of 15 possible triplet motif types were present, as compared to randomized networks. These were in three clear groups: protein-protein interactions, signaling and regulatory triplets where the latter two showed no overlap. In the triplets, negative genetic interactions were associated with paralogs and ohnologs, however these were very rare. Negative genetic interactions among the six triplet motifs did however show strong dosage constraints, with genes being significantly associated with toxicity on overexpression and periodicity in the cell-cycle. Negative genetic interactions overlapped with other interaction types in 37% of cases; these were predominantly associated with protein complexes or signaling events. Finally, we highlight regions of ‘network vulnerability’ containing multiple negative genetic interactions; these could be targeted in fungal species for the regulation of cell growth. 

# Directories


* Negative_Genetic_Interactions
  + Analysis of triplet motifs containing one negative genetic interactions
  
* Random_Edges 
  + Random addition or removal of negative genetic interactions to test the robustness of the six overrepresented triplet motifs
  
* More_stringent_network  
  + We repeated the triplet motif analysis but with data of increasingly stringent genetic interaction scores
  
* Paralogs
  + Pairs of proteins that participate in negative genetic interactions were examined to detect whether they were paralogs or ohnologs

* Analyze_GI_edges
  + Analysis of negative genetic interactions that overlap with some other types of interactions

* Repeated_GI_in_Motifs
  + To find negative genetic interactions that are shared by two or more triplets, which were more frequent than by chance, their frequency was compared between observed and randomized networks

* Negative_GI_Essential
  + To identify triplets that contain significant numbers of essential proteins

* Cell_Cycle
  + To test whether some triplet motifs are enriched for cell cycle-regulated genes 
  
* Overexpressed_Toxic_Genes
  + To test whether some tiplet motifs were enriched for proteins that are toxi upon protein express

* GO_terms
  + To test whether the proteins that share negative genetic interaction are more likely to share the same GO term (GO slim terms)

* Phenotype
  + To test whether the proteins that share negative genetic interaction are more likely to share the same phenotype

* Common
  + Scripts and functions that are commonly shared between multiple scripts

* Figures
  + Scripts that list of all figures (and supplementary figuures) 

* Supplementary_Files
  + Contain souce codes that were commonly used for the analysis
  
* Examples 
  + Some examples for analysis of triplet motifs




* Version 1.0 
* Project start 20th April 2016


* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)