#!/usr/bin/python
import sys
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

dataFile3 = sys.argv[3]
dataSet3 = open(dataFile3)

all_monomers_1 = []
all_monomers_2 = []
MonomerKeys = {}

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

for line in dataSet3.readlines():
    smiles = line.split(" ")[1].strip
    key = line.split(" ")[0]
    # Append to dictionary
    #MonomerKeys.append(key, smiles)

# Generates a list of unique monomers with the number of occurrances of that monomer (Ex: Monomer_A 5)
unique_monos_1 = Counter(all_monomers_1)
unique_monos_1.most_common()

unique_monos_2 = Counter(all_monomers_2)
unique_monos_2.most_common()

# From list of unique monomers, assign a numeric value to each monomer and save in all_monomers_numericKeys list


# Substitue number for SMILES in all_monomers


# Run spearmanr on the two lists of numbers; print (save) spearman values to be reported
#sp.stats.spearmanr([all_monomers_1], [all_monomers_2])


file = open(sys.argv[3], "w")
#for key in all_monomers_1:
#    file.write("%s\n" % key)

for key, count in unique_monos_1.iteritems():
    file.write("%s\t%i\n" % (key, count))

file = open(sys.argv[4], "w")
#for key in all_monomers_2:
#    file.write("%s\n" % key)

for key, count in unique_monos_2.iteritems():
    file.write("%s\t%i\n" % (key, count))

#file.close()