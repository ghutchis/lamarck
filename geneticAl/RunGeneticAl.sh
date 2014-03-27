#!/bin/sh
#$ -S /bin/bash
#$ -R y
#$ -pe mpi 4
#$ -o result.stdout
#$ -e result.stderr
#$ -cwd

#Simple script to submit genetic algorithm to the Sun Grid Engine
#Change the output file name before each run. If the output file name already exists, the text output file won't be saved.
#Seed, length, number in the population (N), matrix nbrs and objective can also be changed as needed


source /etc/profile

export PYTHONPATH=$PYTHONPATH:/nfs/Users/ilanakanal/screeningproject/geneticAl
export PYTHONPATH=$PYTHONPATH:/nfs/Users/ilanakanal/local/lib/python2.7/site-packages

INPUT=sim.inp
WORK=`pwd`
EXE=/nfs/Users/ilanakanal/screeningproject/geneticAl/geneticAl.py
SCRATCH=/scratch/${USER}/screeningproject

# Set up the scratch directory if needed
if [ -d /scratch/${USER} ]; then
 touch /scratch/${USER}
else
 mkdir /scratch/${USER}
fi

mkdir ${SCRATCH}

# Now copy the file into scratch
if [ -f ${INPUT} ]; then
 cp -r ./* ${SCRATCH}
else
exit
fi

cd ${SCRATCH}
python $EXE -d output/March_27_1 --seed=1 --length=6 -N 64 --matrix=mostsim4.json --nbrs=7 --objective=distance


cp -r ${SCRATCH}/* ${WORK}
