#!/bin/bash

#PBS -N TriMoRandEdge[]
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
#PBS -t 1-240


PROJECT_DIRECTORY=Random_Edges
MAIN_SCRIPT=count_triplet_motifs_random_edges.R
HPC_CLUSTER="katana"
NUM_RANDOM_TRIALS=100
LOG_PREFIX="count_triplet_motifs_random_edges_gi"

PERCENTAGE_OF_ORIGINAL_NETWORK=$( seq 0.2 0.1 1.3 )
INNER_LOOP=$(seq 1 1 20)
PARAMS=()

for PERCENTAGE in $PERCENTAGE_OF_ORIGINAL_NETWORK
do

	for OUTPUT in $INNER_LOOP
	do
		echo $PERCENTAGE
		PARAMS+=($PERCENTAGE)
	done
done

# echo ${PARAMS[*]}
#exit

ARRAY_POSITION="$((PBS_ARRAYID - 1))"  

network_size_proportion=${PARAMS[${ARRAY_POSITION}]}


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

if [ ! -d $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${ARRAY_POSITION} ]; then
  mkdir $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${ARRAY_POSITION}
fi

rsync -a ${PBS_O_WORKDIR}/../Common/parameters_file.R  $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${ARRAY_POSITION}
rsync -a ${PBS_O_WORKDIR}/../Common/count_triplet_motifs_helper.R  $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${ARRAY_POSITION}
rsync -a ${PBS_O_WORKDIR}/random_edges_helper.R 			  $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${ARRAY_POSITION}
rsync -a ${PBS_O_WORKDIR}/${MAIN_SCRIPT} 			   $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${ARRAY_POSITION}
rsync -a ${PBS_O_WORKDIR}/../../Data/network_data_library.Rdata    $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${ARRAY_POSITION}

cd $PBS_O_WORKDIR/../../Results/${PROJECT_DIRECTORY}/Job_${ARRAY_POSITION}

echo "My PBS directory"
echo $PBS_O_WORKDIR

echo "My PBS_ARRAYID"
echo ${ARRAY_POSITION}

COMMAND_TO_RUN="Rscript --vanilla ${MAIN_SCRIPT} ${HPC_CLUSTER} ${ARRAY_POSITION} ${NUM_RANDOM_TRIALS} ${network_size_proportion} 1 1 1 > ${LOG_PREFIX}_${ARRAY_POSITION}_${network_size_proportion}_${DATE}.log 2>&1"

echo $COMMAND_TO_RUN
eval $COMMAND_TO_RUN

rm parameters_file.R 
rm count_triplet_motifs_helper.R
rm ${MAIN_SCRIPT}
rm network_data_library.Rdata
rm random_edges_helper.R 


