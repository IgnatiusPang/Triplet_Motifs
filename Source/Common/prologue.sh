#! /bin/bash
 
echo "prologue: copy from working directory to local scratch"
echo "    ${PBS_O_WORKDIR}"
echo "    ${TMPDIR}"
echo
 
#rm -rf ${PBS_O_WORKDIR}/Results
rsync -a ${PBS_O_WORKDIR}/ ${TMPDIR}
