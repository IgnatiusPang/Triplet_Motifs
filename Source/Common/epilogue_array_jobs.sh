#! /bin/bash

echo "epilogue: copy from local scratch to working directory"
echo "    ${TMPDIR}"
echo "    ${PBS_O_WORKDIR}"
echo

rsync -a ${TMPDIR}/*.log ${PBS_O_WORKDIR}/
