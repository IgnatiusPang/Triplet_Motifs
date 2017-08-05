#!/bin/bash

#PBS -N TriMoRandEdge
#PBS -l nodes=1:ppn=16
#PBS -l vmem=16gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M i.pang@unsw.edu.au
#PBS -m ae
#PBS -l epilogue=/srv/scratch/z3371724/Triplet_Motifs/Source/epilogue_random_edges.sh

# Command for running array jobs 
#PBS -t 0-10

### To check he status of array jobs, this is the command:
# qstat  -t -n1 -u z3371724

module load openmpi/1.8.3
module load R/3.2.2

DATE=$(date +"%Y%m%d")

### Create relevant directories

if [ ! -d $PBS_O_WORKDIR/../Results/Random_Edges ]; then
  mkdir $PBS_O_WORKDIR/../Results/Random_Edges
fi

if [ ! -d $PBS_O_WORKDIR/../Log/Random_Edges ]; then
  mkdir $PBS_O_WORKDIR/../Log/Random_Edges
fi

if [ ! -d $PBS_O_WORKDIR/Random_Edges ]; then
  mkdir $PBS_O_WORKDIR/Random_Edges
fi

if [ ! -d $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID} ]; then
  mkdir $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}
fi

PARAMS=( 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.1 1.2 1.3 )

network_size_proportion=${PARAMS[${PBS_ARRAYID}]}

cd $PBS_O_WORKDIR/JDD

rsync -a ${PBS_O_WORKDIR}/count_triplet_motifs_helper.R 	  $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}
rsync -a ${PBS_O_WORKDIR}/random_edges_helper.R 			  $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}
rsync -a ${PBS_O_WORKDIR}/count_triplet_motifs_random_edges.R $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}
rsync -a ${PBS_O_WORKDIR}/network_data_library.Rdata 		  $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}

cd $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}

echo "My PBS directory"
echo $PBS_O_WORKDIR

echo "My PBS_ARRAYID"
echo ${PBS_ARRAYID}

echo "Proportion to Original Size of the Network"
echo $network_size_proportion

Rscript --vanilla count_triplet_motifs_random_edges.R katana ${PBS_ARRAYID} $network_size_proportion 1 1 1 > count_triplet_motifs_random_edges_gi_${network_size_proportion}_$DATE.log

rm $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}/count_triplet_motifs_helper.R 
rm $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}/random_edges_helper.R 
rm $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}/count_triplet_motifs_random_edges.R 
rm $PBS_O_WORKDIR/Random_Edges/Job_${PBS_ARRAYID}/network_data_library.Rdata 
