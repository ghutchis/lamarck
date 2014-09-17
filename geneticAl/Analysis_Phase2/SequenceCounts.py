#!/usr/bin/python
import sys
from string import Template
from collections import Counter
import numpy as np
import scipy as sp

# Manually convert .db file to text file. Data files 1,2 are the output databases which have been converted to text files.
# Data file 3 is the output file name.

dataFile1 = sys.argv[1]
dataSet1 = open(dataFile1)

dataFile2 = sys.argv[2]
dataSet2 = open(dataFile2)

all_sequences_1 = []
all_sequences_2 = []
#TODO make sure to retain the generation throughout this process so that analysis can be done to see how many generations needed for convergence
for line in dataSet1.readlines():
    # pull the seuqnces
    sequence = line.split(",")[2].split("\"")[1].split("\"")[0]
    all_sequences_1.append(sequence)

for line in dataSet2.readlines():
    # pull the two smiles
    smiles = line.split("\"")[1].split("_")[0].split("~")
    all_monomers_2.append(smiles[0])
    all_monomers_2.append(smiles[1])


all_monomers_1.sort()
all_monomers_2.sort()

# Generates a list of unique monomers with the number of occurrances of that monomer (Ex: Monomer_A 5)
unique_monos_1 = Counter(all_monomers_1)
unique_monos_1.most_common()

unique_monos_2 = Counter(all_monomers_2)
unique_monos_2.most_common()


file = open(sys.argv[4], "w")
#for key in all_monomers_1:
#    file.write("%s\n" % key)

run = 1

for key, count in unique_monos_1.iteritems():
    file.write("%i\t%s\t%i\n" % (run, key, count))

#for key in all_monomers_2:
#    file.write("%s\n" % key)

run = 2

for key, count in unique_monos_2.iteritems():
    file.write("%i\t%s\t%i\n" % (run, key, count))

file.close()