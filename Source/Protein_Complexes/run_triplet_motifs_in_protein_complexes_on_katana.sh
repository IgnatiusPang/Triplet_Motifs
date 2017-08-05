#!/bin/bash

#PBS -N TriMoComplex
#PBS -l nodes=1:ppn=16
#PBS -l vmem=8gb
#PBS -l walltime=11:00:00
#PBS -j oe
#PBS -M i.pang@unsw.edu.au
#PBS -m ae
#PBS -l epilogue=/srv/scratch/z3371724/Triplet_Motifs/Source/epilogue.sh

module load openmpi/1.8.3
module load R/3.2.2

DATE=$(date +"%Y%m%d")

cd $PBS_O_WORKDIR

rsync -a ${PBS_O_WORKDIR}/count_triplet_motifs_helper.R ${TMPDIR}
rsync -a ${PBS_O_WORKDIR}/count_triplet_motifs_protein_complexes.R ${TMPDIR}
rsync -a ${PBS_O_WORKDIR}/network_data_library.Rdata ${TMPDIR}

cd $TMPDIR

Rscript --vanilla  count_triplet_motifs_protein_complexes.R katana > count_triplet_motif_protein_complexes_$DATE.log


