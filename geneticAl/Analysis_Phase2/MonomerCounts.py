#!/usr/bin/python
import sys
from string import Template
from collections import Counter
import numpy as np
import scipy as sp

# Manually convert .db file to text file. Data files 1,2 are the output databases which have been converted to text files.
# Data file 3 is the <dims_and_tets folder>_mono.txt file which is used for SMILES substitution. Move this file from the
# dims_and_tets folder or make sure to call it from that location.
# Data file 4 is the output file name.

dataFile1 = sys.argv[1]
dataSet1 = open(dataFile1)

dataFile2 = sys.argv[2]
dataSet2 = open(dataFile2)

all_monomers_1 = []
all_monomers_2 = []

for line in dataSet1.readlines():
    # pull the two smiles
    smiles = line.split("\"")[1].split("_")[0].split("~")
    all_monomers_1.append(smiles[0])
    all_monomers_1.append(smiles[1])

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

# Save the monomer frequencies to a file as defined when running the program
file = open(sys.argv[3], "w")

run = 1

for key, count in unique_monos_1.iteritems():
    file.write("%i\t%s\t%i\n" % (run, key, count))

run = 2

for key, count in unique_monos_2.iteritems():
    file.write("%i\t%s\t%i\n" % (run, key, count))

file.close()