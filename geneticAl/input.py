#!/usr/bin/python

import sys
import pybel
import openbabel

def globalopt(mol, debug=False, fast=False):
    pybel._builder.Build(mol.OBMol)
    mol.addh()

    ff = pybel._forcefields["mmff94"]
    success = ff.Setup(mol.OBMol)
    if not success:
        ff = pybel._forcefields["uff"]
        success = ff.Setup(mol.OBMol)
        if not success:
            sys.exit("Cannot set up forcefield")

    ff.SteepestDescent(1000, 1.0e-4)    
    ff.WeightedRotorSearch(250, 25)
    ff.WeightedRotorSearch(250, 25)
    ff.ConjugateGradients(500, 1.0e-6)
    ff.GetCoordinates(mol.OBMol)

header = "%%nproc=4\n#T PM6 OPT(MaxCycle=1000)"
header_b = "\n\n"

datafile = open(sys.argv[1])

for line in datafile:
    line = line.rstrip()
    name = line.split()[0]
    #name = name.replace(':', '-')
    dimer = line.split()[0] * 2
    #tetramer = line.split()[0] * 4
    #mol = pybel.readstring("smi", tetramer)
    mol = pybel.readstring("smi", dimer)
    print name

    globalopt(mol)

    gaussian = (header + "\n\n" + name + "\n" + "\n".join(mol.write("gau").replace("0  3\n", "0  1\n").split("\n")[3:]))

    output = open("%s.com" % name, "w")
    output.write(gaussian)
    output.close()