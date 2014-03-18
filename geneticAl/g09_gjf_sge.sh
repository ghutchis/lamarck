#!/bin/sh
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
  cp ${INPUT%%.gjf}.* ${SCRATCH}
else
  exit
fi

# Set up G09 environment
export g09root=/usr/local
export GAUSS_SRCDIR=${SCRATCH}
. $g09root/g09/bsd/g09.profile


cd ${SCRATCH}

# If present, re-generate the chk file
fchk=${INPUT%%.gjf}.fchk
if [ -f ${fchk}.gz ]; then
  gunzip ${fchk}.gz
fi
if [ -f ${fchk} ]; then
  unfchk ${fchk}
fi

# Remove any old *.gz files
if [ `pwd` -ef ${SCRATCH} ]; then
  rm *.gz
fi
FILENAME=${INPUT%%.gjf}
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