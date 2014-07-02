#!/usr/bin/python
import sys
from collections import Counter

# Convert .db file to text file. Open in excel and save the first column as a new document in .csv format
dataFile = sys.argv[1]
dataSet = open(dataFile)

all_monomers = []
unique_moomers = []

for line in dataSet.readlines():
    # pull the two smiles
    smiles = line.split("\"")[1].split("_")[0].split("~")
    print smiles
    all_monomers.append(smiles[0])
    all_monomers.append(smiles[1])

all_monomers.sort()
print all_monomers