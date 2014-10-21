#!/usr/bin/python
import sys
import json
import numpy as np
import pandas as pd

# load dataset
dataFile1 = sys.argv[1]
dataSet1 = open(dataFile1)

# expand text string to identify SMILES, sequence, and other stuff
monomers = {}

for line in dataSet1.readlines():
    # pull the two smiles
    smiles = line.split("\"")[1].split("_")[0].split("~")
    sequence = line.split(",")[2].split("\"")[1].split("\"")[0]

    # create dict of dict based on monomer & sequence, with values being counts
    # if the monomer and sequence exists...
    for mon in smiles:
        try:
            # if both monomer and sequence exist, increment by 1
            monomers[mon][sequence] += 1
        except:
            try:
                # if sequence doesn't exists but monomer does, increment by 1
                monomers[mon][sequence] = 1
            except:
                # if neither exist, create
                monomers[mon] = {sequence: 1}

# save the dictdict (see http://stackoverflow.com/a/7100202/168775)
with open('SequenceCounts.json', 'wb') as fp:
    json.dump(monomers, fp)

# replace the NaNs with 0 (see http://stackoverflow.com/a/13295801/168775)
monomerDF = pd.DataFrame(monomers).transpose().fillna(0)

monomerDF.to_csv('SequenceCounts.csv')