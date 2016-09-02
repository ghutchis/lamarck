#!/usr/bin/python
import sys
from string import Template
from collections import Counter
import numpy as np
import scipy as sp
import pybel
import openbabel
import csv

# Manually convert .db file to text file. Data file is a file which contains the output databases which have been
# converted to text files and merged into one large file.

dataFile1 = sys.argv[1]
dataSet1 = open(dataFile1)

all_monomer_pairs = []

for line in dataSet1.readlines():
    # Pull SMILES
    smiles = line.split("_")[0].split("~")
    #smiles = line.split("\"")[1].split("_")[0].split("~")
    # Convert SMILES to canonical SMILES
    mol_1 = pybel.readstring("smi", smiles[0])
    canmol1 = mol_1.write("can").split("\t")[0]
    mol_2 = pybel.readstring("smi", smiles[1])
    canmol2 = mol_2.write("can").split("\t")[0]
    # Make a set containing the 2 canonical SMILES and save the set to a list of all monomer pairs
    all_monomer_pairs.append({canmol1, canmol2})

# Generates a list of unique pairs with the number of occurrances of that monomer (Ex: Monomer_A 5)
unique_data = [list(x) for x in set(frozenset(tuple(x)) for x in all_monomer_pairs)]

monomer_pair_counts = []
for pair in unique_data:
    counts = all_monomer_pairs.count(set(pair))
    monomer_pair_counts.append([pair,counts])

sorted_pairs = sorted(monomer_pair_counts, key=lambda tup: tup[1], reverse=True)
# Save the monomer frequencies to a file as defined when running the program
with open(sys.argv[2], "wb") as f:
    writer = csv.writer(f)
    writer.writerows(sorted_pairs)
