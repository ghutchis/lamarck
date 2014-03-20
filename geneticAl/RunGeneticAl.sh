#!/bin/sh

#Simple script to submit genetic algorithm to the Sun Grid Engine
#Change the output file name before each run. If the output file name already exists, the text output file won't be saved.
#Seed, length, number in the population (N), matrix nbrs and objective can also be changed as needed

# Force SGE to use bash
#$ -S /bin/bash
#$ -R y
#$ -pe mpi 4

# Take the first argument of the script
INPUT=$1
BASE=$2
SCRATCH=/scratch/${USER}/${JOB_ID}

# Set up the scratch directory if needed
if [ -d /scratch/${USER} ]; then
  touch /scratch/${USER}
else
  mkdir /scratch/${USER}
fi
mkdir ${SCRATCH}

# Now copy the file into scratch
cd $BASE
if [ -f ${INPUT} ]; then
  # also copy any fchk or chk file with the same name
  cp ${INPUT%%.py}.* ${SCRATCH}
  cp ${INPUT%%.gjf}.* ${SCRATCH}
else
  exit
fi


qsub -N python  geneticAl.py -d output/Run3_18_1 --seed=1 --length=6 -N 6 --matrix=mostsim4.json --nbrs=7 --objective=distance