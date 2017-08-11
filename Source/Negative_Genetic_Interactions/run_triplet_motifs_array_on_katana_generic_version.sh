#!/bin/bash

#PBS -N TriMotif[]
#PBS -l nodes=1:ppn=16
#PBS -l vmem=100gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M i.pang@unsw.edu.au
#PBS -m ae
## Clive
#### PBS -l epilogue=/srv/scratch/z3371724/Triplet_Motifs/Source/epilogue_array_jobs.sh

## Katana
#PBS -l epilogue=/srv/scratch/z3371724/Triplet_Motifs/Source/epilogue_array_jobs.sh

# Command for running array jobs 
#PBS -t 1-20

PROJECT_DIRECTORY=Negative_Genetic_Interactions
MAIN_SCRIPT=count_triplet_motifs_negative_interactions.R
HPC_CLUSTER="clive"
NUM_RANDOM_TRIALS=100
LOG_PREFIX="count_triplet_motifs_neg_int"

if [ $HPC_CLUSTER = 'clive' ] 
then 
	## Clive
	module load openmpi/1.10.3
	module load R/3.3.1-bioconductor-3.3
elif [ $HPC_CLUSTER = 'katana' ]
then 
	## Katana
	module load openmpi/1.8.3
	module load R/3.2.2
else
	echo "Invalid HPC_CLUSTER = ${HPC_CLUSTER}"
fi


### To check he status of array jobs, this is the command:
# qstat  -t -n1 -u z3371724

DATE=$(date +"%Y%m%d")


export PROJECT_DIRECTORY

### Create relevant directories

if [ ! -d $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY} ]; then
  mkdir $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}
fi

if [ ! -d $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${PBS_ARRAYID} ]; then
  mkdir $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${PBS_ARRAYID}
fi

rsync -a ${PBS_O_WORKDIR}/../Common/parameters_file.R  $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${PBS_ARRAYID}
rsync -a ${PBS_O_WORKDIR}/../Common/count_triplet_motifs_helper.R  $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${PBS_ARRAYID}
rsync -a ${PBS_O_WORKDIR}/${MAIN_SCRIPT} 			   $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${PBS_ARRAYID}
rsync -a ${PBS_O_WORKDIR}/../../Data/network_data_library.Rdata    $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${PBS_ARRAYID}

cd $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${PBS_ARRAYID}

echo "My PBS directory"
echo $PBS_O_WORKDIR

echo "My PBS_ARRAYID"
echo ${PBS_ARRAYID}

COMMAND_TO_RUN="Rscript --vanilla ${MAIN_SCRIPT} ${HPC_CLUSTER} ${PBS_ARRAYID} ${NUM_RANDOM_TRIALS} > ${LOG_PREFIX}_${PBS_ARRAYID}_${DATE}.log 2>&1"

echo $COMMAND_TO_RUN

`$COMMAND_TO_RUN`

rm parameters_file.R 
rm count_triplet_motifs_helper.R
rm ${MAIN_SCRIPT}
rm network_data_library.Rdata


