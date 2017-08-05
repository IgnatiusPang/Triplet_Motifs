
## Transfer scripts on testing negative genetic interactions from Rm244 to clusters (array job )


### Clive

rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Data/Triplet_Motifs_R_data/network_data_library.Rdata z3371724@clive:/share/bioinfo/igy/Triplet_Motifs/Data


## Send the common scripts over to clive
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Common/count_triplet_motifs_helper.R z3371724@clive:/share/bioinfo/igy/Triplet_Motifs/Source
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Common/parameters_file.R             z3371724@clive:/share/bioinfo/igy/Triplet_Motifs/Source


rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Negative_Genetic_Interactions/count_triplet_motifs_negative_interactions.R    z3371724@clive:/share/bioinfo/igy/Triplet_Motifs/Source/Negative_Genetic_Interactions/
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Negative_Genetic_Interactions/run_triplet_motifs_array_on_clive.sh z3371724@clive:/share/bioinfo/igy/Triplet_Motifs/Source/Negative_Genetic_Interactions/
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Negative_Genetic_Interactions/epilogue_array_jobs.sh z3371724@clive:/share/bioinfo/igy/Triplet_Motifs/Source/Negative_Genetic_Interactions/

## Testing scripts
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Common/parameters_file_small_tests.R             z3371724@clive:/share/bioinfo/igy/Triplet_Motifs/Source/Negative_Genetic_Interactions/
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Negative_Genetic_Interactions/run_triplet_motifs_array_on_clive_generic_version.sh z3371724@clive:/share/bioinfo/igy/Triplet_Motifs/Source/Negative_Genetic_Interactions/
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Negative_Genetic_Interactions/epilogue_array_jobs_testing.sh z3371724@clive:/share/bioinfo/igy/Triplet_Motifs/Source/Negative_Genetic_Interactions/


### Katana
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Negative_Genetic_Interactions/count_triplet_motifs_negative_interactions.R   z3371724@kdm.science.unsw.edu.au:/srv/scratch/z3371724/Triplet_Motifs/Source/
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Negative_Genetic_Interactions/run_triplet_motifs_array_on_katana.sh z3371724@kdm.science.unsw.edu.au:/srv/scratch/z3371724/Triplet_Motifs/Source/
rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Source/Negative_Genetic_Interactions/epilogue_array_jobs.sh z3371724@kdm.science.unsw.edu.au:/srv/scratch/z3371724/Triplet_Motifs/Source/




## Get results

rsync -avh /media/z3371724/PostDoc/2016/Triplet_Motifs/Data/Triplet_Motifs_R_data/network_data_library.Rdata z3371724@kdm.science.unsw.edu.au:/srv/scratch/z3371724/Triplet_Motifs/Data


## From Katana
rsync -avh  z3371724@kdm.science.unsw.edu.au:/srv/scratch/z3371724/Triplet_Motifs/Source/Neg_Genetic_Interactions/Job_*  /media/z3371724/PostDoc/2016/Triplet_Motifs/Results/Bootstrap_p_values/Negative_Genetic_Interactions/


