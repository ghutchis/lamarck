__author__ = 'ilanakanal'
#!/usr/bin/python
import sys
import numpy as np
import subprocess



# get number of lines in a file
def file_len(f):
    f = open(f)
    for i, l in enumerate(f):
        pass
    return i + 1


# find number of generations in a results file
def get_num_generations(f):
    l = subprocess.check_output(['tail', '-1', f])
    l = l.split(",")[1]
    return l[1:-1]

# todo: eventually, we want to include a loop here to loop over any number of input files provided to the program and
# todo: output the whole thing as a multidim array. Or a number of files.
# todo: Whatever.

# datafile format is .db from geneticAl.py saved as a .txt file with the header line removed
dataFile = sys.argv[1]
monomerFile = "/Users/ilanakanal/Code/screeningproject/monomers-300b.txt"

rows = int(file_len(monomerFile))
cols = int(get_num_generations(dataFile)) + 1  # generations start at 0
monomerCounts = np.zeros((rows, cols))

dataSet = open(dataFile)
monomerSet = open(monomerFile).read().splitlines()

for n, smile in monomerSet:
    monomerSet(n) = smile.strip()

print monomerSet
for line in dataSet.readlines():
    # pull the two smiles
    smiles = line.split("_")[0][1:].split("~")
    # pull the generation
    gen = int(line.split(",")[1][1:-1])

    for monomer in smiles:
        # find smile location in monomerSet
        monomerCounts[monomerSet.index(monomer), gen] += 1

# save it out
np.savetxt("stuff.csv", monomerCounts, delimiter=",")

print("""
The output is saved as stuff.csv. THERE ARE NO ROW LABELS; simply add a column at the beginning
with the contents of monomerfile-300b.txt. The output data is unsorted. Sum across the rows and
sort descending to see most common monomers.
""")