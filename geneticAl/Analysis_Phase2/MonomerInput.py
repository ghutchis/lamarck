#!/usr/bin/python

# Makes monomer gaussian input files to use in PCA analysis for geneticAl. Before beginning, make a directory in which
# to save the .com files. For running this file, sys.argv[1] is the monomer file saved in the folder in dims_and_tets
# from generating the similarity matrix.

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

    ff.SteepestDescent(500, 1.0e-4)
    ff.WeightedRotorSearch(200, 25)
    ff.ConjugateGradients(500, 1.0e-6)
    ff.GetCoordinates(mol.OBMol)

header = "%%nproc=4\n%%Chk=%s.chk\n#T B3LYP/6-31G* OPT(MaxCycle=1000)"
header_b = "\n\n"

datafile = open(sys.argv[1])


for line in datafile:
    line = line.rstrip()
    name = line.split()[0]
    mol = pybel.readstring("smi", line.split()[1])
    print name

    globalopt(mol)

    gaussian = (header + "\n\n" + name + "\n" +
                "\n".join(mol.write("gau").replace("0  3\n", "0  1\n").split("\n")[3:])) % (name)

    output = open("%s.com" % name, "w")
    output.write(gaussian)
    output.close()
