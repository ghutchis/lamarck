#!/usr/bin/python
import sys
from collections import Counter


# Convert .db file to text file.
dataFile = sys.argv[1]
dataSet = open(dataFile)

all_monomers = []

for line in dataSet.readlines():
    # pull the two smiles
    smiles = line.split("\"")[1].split("_")[0].split("~")
    all_monomers.append(smiles[0])
    all_monomers.append(smiles[1])

all_monomers.sort()

unique_monos = Counter(all_monomers)
unique_monos.most_common()

file = open(sys.argv[2], "w")

for key, count in unique_monos.iteritems():
    file.write("%s\t%i\n" % (key, count))

file.close()