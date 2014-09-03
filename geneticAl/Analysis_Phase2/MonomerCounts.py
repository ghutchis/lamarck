#!/usr/bin/python
import sys
from collections import Counter


# Manually convert .db file to text file.

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

unique_monos_1 = Counter(all_monomers_1)
unique_monos_1.most_common()

unique_monos_2 = Counter(all_monomers_1)
unique_monos_2.most_common()


#file = open(sys.argv[3], "w")

#for key, count in unique_monos.iteritems():
#    file.write("%s\t%i\n" % (key, count))

#file.close()