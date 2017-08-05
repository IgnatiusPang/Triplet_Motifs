#! /bin/bash

echo "epilogue: copy from local scratch to working directory"
echo "    ${TMPDIR}"
echo "    ${PBS_O_WORKDIR}"
echo

rsync -a ${TMPDIR}/*.tab ${PBS_O_WORKDIR}/../Results

rsync -a ${TMPDIR}/*.Rdata ${PBS_O_WORKDIR}/../Results

rsync -a ${TMPDIR}/*.log  ${PBS_O_WORKDIR}/../Log


