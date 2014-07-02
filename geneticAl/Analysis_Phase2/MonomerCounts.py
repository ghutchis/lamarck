#!/usr/bin/python
import sys
from collections import Counter

# Convert .db file to text file. Open in excel and save the first column as a new document in .csv format
dataFile = sys.argv[1]
dataSet = open(dataFile)

for line in dataSet.readlines():
    # pull the two smiles
    smiles = line.split("_")[0].split("~")
    print smiles

    #for item in smiles:
     #   str = "-"
      #  smilesList = str.join(smiles)
       # print smilesList

  #  for monomer in smiles:
