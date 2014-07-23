import os
import sys
import json
import math

def getdata(folder):
    """Read the zindo.txt spreadsheet for the homopolymers"""
    input = open(os.path.join(folder, "zindo.txt"))
    header = input.next()
    data = []
    for line in input:
        broken = line.split()
        file, id, smile, homo, lumo = broken[:5]
        homo, lumo = float(homo), float(lumo)
        trans = map(float, broken[5:])
        etenergy = []
        osc = []
        for i in range(0, len(trans), 2):
            etenergy.append(trans[i])
            osc.append(trans[i+1])
        data.append((smile, homo, lumo, etenergy, osc))
    return data

def makesimmatrix(folder, length):
    SMILE, HOMO, LUMO, ETENERGY, OSC = range(5)
    lumos = []
    etenergies = []

    data = getdata(folder)
    cutoff = 0.1
    newdata = {}
    for d in data:
        chosen = min([x for (x,y) in enumerate(d[OSC]) if y>=cutoff] + [999])
        if chosen == 999:
            mymax = max(d[OSC])
            chosen = [x for (x,y) in enumerate(d[OSC]) if y==mymax][0]
        smile = d[SMILE][:(len(d[SMILE]) / length)]
        newdata[smile] = (d[LUMO], d[LUMO] - d[HOMO])
        
    smiles = newdata.keys()
    smiles.sort()
    lookup = dict((y,x) for x,y in enumerate(smiles))
    simmat = [[0]*len(smiles) for x in smiles]

    mostsimilar = {}
    for smile in smiles:
        i = lookup[smile]
        d = newdata[smile]
        for smile_b in smiles:
            j = lookup[smile_b]
            e = newdata[smile_b]
            simmat[i][j] = math.sqrt((d[0] - e[0])**2 + (d[1]-e[1])**2)
        mostsimilar[smile] = [y for x,y in sorted(zip(simmat[i], smiles))
                              if y!=smile]
    
    debug = False
    if debug:
        print "(05/04/2012) This looks like debugging output to me, to check"
        print "that the code is selecting the mostsimilar monomer in each case."
        for idx in range(0, 2):
            print "\n\n\n", idx
            print mostsimilar[smiles[idx]]  
            minidx = [x for x,y in enumerate(simmat[idx]) if y==min(simmat[idx][(idx+1):])][0]
            print minidx, simmat[idx][minidx]
            print smiles[idx]
            print newdata[smiles[idx]]
            print smiles[minidx]
            print newdata[smiles[minidx]]


    # What is the smallest number of similar components
    # that includes all molecules
    found = -1
    for i in range(5, 40):
        allval = set()
        for v in mostsimilar.values():
            allval.update(v[:i])
        if len(mostsimilar) == len(allval) and found == -1:
            found = i
##        print i, len(allval), len(mostsimilar), set(mostsimilar.keys()) - allval
    assert found != -1
    
    return mostsimilar, found

def usage():
    print """
  python MakeSimMatrix.py folder length outputfile
  
    where
      folder should contain a zindo.txt from running ExtractData.py
      length is the length of the polymer in repeat units
        (this value should be the same as you used when running
         Homopolymer.py)
      outputfile will be a JSON file containing the list of monomers
        and their associated neighbours (in order of similarity).
        
  """
    
    sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()

    folder = sys.argv[1]
    if not os.path.isdir(folder) or not os.path.isfile(os.path.join(folder, "zindo.txt")):
        usage()
    length = int(sys.argv[2])
    outputfilename = sys.argv[3]

    mostsim, found = makesimmatrix(folder, length)
    with open(outputfilename, "w") as f:
        json.dump(mostsim, f)

    print """

***************************************************************
When running the Genetic Algorithm using this similarity matrix,
you should use at least
     %d   nearest neighbours (MAKE A NOTE!!)
or else some of your monomers will never be mutated to.

Creating mostsim.json
***************************************************************""" % found

##monolen2 136 8 27
##monolen4 133 50 63
##monolen6 134 99 111
##monolen8 132 106 117
    
