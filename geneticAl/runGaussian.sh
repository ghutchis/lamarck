#!/bin/sh
#
# Simple script to submit gjf files from genetic algorithm to the Sun Grid Engine

# Take the first argument of the script
BASE=`pwd`

for x in "$@"; do
 INPUT=${x}
 qsub -N ${INPUT%%.gjf} /usr/local/bin/g09_gjf_sge.sh ${INPUT} ${BASE}
 echo "Gaussian job ${INPUT} submitted to the queue."
done