#!/bin/sh

#Simple script to submit genetic algorithm to the Sun Grid Engine
#Change the output file name before each run. If the output file name already exists, the text output file won't be saved.
#Seed, length, number in the population (N), matrix nbrs and objective can also be changed as needed

python  geneticAl.py -d output/Run3_9_1 --seed=1 --length=6 -N 6 --matrix=mostsim4.json --nbrs=7 --objective=distance