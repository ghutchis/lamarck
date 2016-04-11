#!/bin/sh

#PBS -N GeneticAl
#PBS -j oe
#PBS -l nodes=1:ppn=12
#PBS -q dev
#PBS -l walltime=00:14:00

# Simple script to submit genetic algorithm to Frank
# Output file name should be specified when call this to run. If the output file name already exists, the text output
# file won't be saved.
# Seed, length, number in the population (N), matrix nbrs and objective can also be changed as needed

# Each new user to run the genetic algortihm on Frank must install in their own user space cclib as described below,
# (make sure you have the desired python loaded, I tested with python/2.7-gcc45):
# mkdir python-modules
# git clone https://github.com/cclib/cclib.git
# cd cclib
# python setup.py install --user

cd $SCRATCH

module load python/2.7-gcc45
module load openbabel
module load cclib/1.1
module load gaussian/g09B.01

# INPUT=$1
WORK=`pwd`
EXE=/ihome/ghutchison/iek2/screeningproject/geneticAl/geneticAl.py
#SCRATCH=/scratch/${USER}/${JOB_ID}

python $EXE -d output/Apr_11_16_a --seed=2 --length=4 -N 6 --matrix=mostsim4.json --nbrs=7 --objective=distance

cp $SCRATCH/* $PBS_O_WORKDIR
