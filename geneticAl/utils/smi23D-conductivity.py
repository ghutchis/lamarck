"""
Calculate conductivity
"""
import os
import sys
import pdb

import pybel
try: # Works for Python 2.6+
    import json as simplejson
except ImportError: # Otherwise you need to install simplejson
    import simplejson

from Utils import molname_to_mol, globalopt

def exitwitherror():
    sys.exit(__doc__)

header = "%%nproc=4\n%%mem=1GB\n%%Chk=%s.chk\n#T B3LYP/6-31G(d) OPT"
header_cation = """
%%nproc=4
%%mem=1GB
%%Chk=%s.chk
%%NoSave
#T B3LYP/6-31G(d) Opt Geom=Checkpoint Guess=Read

Cation %s

1  2

"""
header_singlepoint = "%nproc=1\n%mem=1GB\n#T B3LYP/6-31G(d) SP"

if __name__ == "__main__":
    if len(sys.argv) < 3:
        exitwitherror()

    idx = sys.argv[1]
    length = int(sys.argv[2])

    input = open(os.path.join("gaussian", "%s.txt" % idx), "r")
    smi = input.read().rstrip()
    input.close()

    mol = molname_to_mol(smi, length)
    mol.make3D()
    globalopt(mol)

    # basic idea:
    # B3LYP geometry optimization of neutral, cation (above)
    # Single point energy of neutral@cation geometry, cation@neutral geometry
    # Take all 4 SCF energies
    # "Reorganization energy" = (neutral@cation - neutral) + (cation@neutral geometry - cation)

    gaussian = (header % idx + "\n\n" + smi + "\n" + 
                "\n".join(mol.write("gau").replace("0  3\n", "0  1\n").split("\n")[3:]))
    inputFile = os.path.join("gaussian", "%s.gjf" % idx)
    outputFile = os.path.join("gaussian", "%s.out" % idx)

    output = open(inputFile, "w")
    output.write(gaussian)
    output.close()

    # Run gaussian calculation to get neutral geometry
    os.system("/usr/local/g09/g09 %s %s" %(inputFile, outputFile))

    file = pybel.readfile("g03", outputFile)
    neutralMol = file.next()
    neutralGeomE = neutralMol.energy

    # OK, neutral is done, now run cation optimization
    gaussian = (header_cation % (idx, smi))
    inputFile = os.path.join("gaussian", "%s+.gjf" % idx)
    outputFile = os.path.join("gaussian", "%s+.out" % idx)

    output = open(inputFile, "w")
    output.write(gaussian)
    output.close()

    # Uses checkpoint of neutral to save time
    os.system("/usr/local/g09/g09 %s %s" %(inputFile, outputFile))
    
    file = pybel.readfile("g03", outputFile)
    cationMol = file.next()
    cationGeomE = cationMol.energy

    # Run single point of neutral@cation geometry
    gaussian = (header_singlepoint + "\n\n SP of Neutral @ Cation geometry" + smi + "\n" + 
                "\n".join(cationMol.write("gau").replace("1  2\n", "0  1\n").split("\n")[3:]))
    inputFile = os.path.join("gaussian", "%s@.gjf" % idx)
    outputFile = os.path.join("gaussian", "%s@.out" % idx)

    output = open(inputFile, "w")
    output.write(gaussian)
    output.close()

    os.system("/usr/local/g09/g09 %s %s" %(inputFile, outputFile))
    file = pybel.readfile("g03", outputFile)
    neutralAtCat = file.next()
    neutralAtCationE = neutralAtCat.energy

    # Run single point of cation@neutral geometry
    gaussian = (header_singlepoint + "\n\n SP of Cation @ Neutral geometry" + smi + "\n" + 
                "\n".join(neutralMol.write("gau").replace("0  1\n", "1  2\n").split("\n")[3:]))
    inputFile = os.path.join("gaussian", "%s+@.gjf" % idx)
    outputFile = os.path.join("gaussian", "%s+@.out" % idx)

    output = open(inputFile, "w")
    output.write(gaussian)
    output.close()

    os.system("/usr/local/g09/g09 %s %s" %(inputFile, outputFile))
    file = pybel.readfile("g03", outputFile)
    catAtNeutral = file.next()
    catAtNeutralE = catAtNeutral.energy

    print smi, length, (neutralAtCationE - neutralGeomE) + (catAtNeutralE - cationGeomE)
