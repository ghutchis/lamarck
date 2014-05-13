#!/bin/sh
#$ -S /bin/bash
#$ -R y
#$ -pe mpi 4
#$ -o result.stdout
#$ -e result.stderr
#$ -cwd

# Simple script to submit genetic algorithm to the Sun Grid Engine
# Output file name should be specified when call this to run. If the output file name already exists, the text output
# file won't be saved.
# Seed, length, number in the population (N), matrix nbrs and objective can also be changed as needed

source /etc/profile

#set up python and g09 environments
export PYTHONPATH=$PYTHONPATH:/nfs/Users/ilanakanal/screeningproject/geneticAl
export PYTHONPATH=$PYTHONPATH:/nfs/Users/ilanakanal/local/lib/python2.7/site-packages
#export g09root=/usr/local
#export GAUSS_SRCDIR=${SCRATCH}
. $g09root/g09/bsd/g09.profile

# INPUT=$1
WORK=`pwd`
EXE=/nfs/Users/ilanakanal/screeningproject/geneticAl/geneticAl.py
SCRATCH=/scratch/${USER}/${JOB_ID}

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
python $EXE -d output/March_31_2 --seed=2 --length=4 -N 6 --matrix=mostsim4.json --nbrs=7 --objective=distance


cp -r ${SCRATCH}/* ${WORK}

if [ `pwd` -ef ${SCRATCH} ]; then
  rm *.gz
fi
FILENAME=${INPUT%%.com}
g09 ${INPUT} ${FILENAME}.g09

# Cleanup
if [ -f ${FILENAME}.g09 ]; then
  gzip -9 ${FILENAME}.g09
  cp ${FILENAME}.g09.gz ${BASE}
fi

if [ -f ${FILENAME}.chk ]; then
  formchk -3 ${FILENAME}.chk
elif [ `echo *.chk` != '*.chk' ]; then
  mv `echo *.chk` ${FILENAME}.chk
  formchk -3 ${FILENAME}.chk
fi

if [ -f ${FILENAME}.fchk ]; then
  gzip -9 ${FILENAME}.fchk
  cp ${FILENAME}.fchk.gz ${BASE}
fi

# Don't delete scratch directory for now
cd $BASE
rm -rf ${SCRATCH}