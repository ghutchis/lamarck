# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:06:31 2014

@author: ilanakanal
"""

"""
"""
import sys
import pdb
import os
import pybel
import json

#from Utils import molname_to_mol, globalopt
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

def _readdatafiles():
    relpath = os.sep.join(__file__.split(os.sep)[:-1])

    f = open(os.path.join(relpath, "directions.json"))
    # Need to remove the unicode from the json or else Open Babel won't like it
    combinations = dict([(float(x), [[str(t) for t in z] for z in y]) for x,y in json.load(f).iteritems()])
    f.close()
    return combinations
    
combinations = _readdatafiles()

def getreversed(smiles):
    mol = pybel.readstring("smi", smiles)
    dimer = pybel.readstring("smi", smiles+smiles)

    # Identify the 'last' atom (not necessarily the last
    # atom in the SMILES string)
    # One way to do this is to find what atom has a new
    # connection in the dimer
    last = 0
    for i, (a,b) in enumerate(zip(mol, dimer)):
        if a.OBAtom.GetValence() != b.OBAtom.GetValence():
            last = i+1
            break
    assert last != 0

    # Generate the reverse
    conv = pybel.ob.OBConversion()
    conv.SetOutFormat("smi")
    conv.SetOptions('f"%d"l"%d"' % (last, 1), conv.OUTOPTIONS)
    reverse = conv.WriteString(mol.OBMol).rstrip()

    # Sanity check!! Dimer of 'reverse' should be identical
    newdimer = pybel.readstring("smi", reverse+reverse)
    a, b = dimer.write("can"), newdimer.write("can")
    assert a == b, "\n %s is not the same as\n %s" % (a, b)
    
    return reverse
    
def issym(smiles):
    """Is this monomer symmetric?
    """
    rev = getreversed(smiles)
    return (pybel.readstring("smi", smiles+"Br").write("can") ==
            pybel.readstring("smi", rev + "Br").write("can"))

def molname_to_smile(molname, length):
    """Convert a molecule 'name' to a molecule

    >>> molname = "c(s1)c(OCO2)c2c1~C=CN=N_1_0"
    >>> smile = molname_to_smile(molname, 4)
    >>> print smile
    C=CN=Nc(s1)c(OCO2)c2c1C=CN=Nc(s1)c(OCO2)c2c1
    """
    parts, d1, d2 = molname.split("_")
    b = parts.split("~") # monos

    comb = combinations[length][int(d1)][int(d2)]
    reversed = []
    for i in range(2):
        if issym(b[i]):
            reversed.append(b[i])
        else:
            reversed.append(getreversed(b[i]))
    smile = comb.replace("(qu)", b[0]).replace("(uq)", reversed[0]).replace(
                         "(di)", b[1]).replace("(id)", reversed[1])
    return smile

def molname_to_mol(molname, length):
    smi = molname_to_smile(molname, length)
    return pybel.readstring("smi", smi)

def exitwitherror():
    sys.exit(__doc__)

header = "%%nproc=4\n%%mem=1GB\n%%Chk=%s.chk\n#T PM6 OPT(MaxCycle=1000)"
header_b = """
--Link1--
%%nproc=4
%%mem=1GB
%%Chk=%s.chk
%%NoSave
# Geom=AllCheck ZINDO(NStates=15,Singlets)
"""

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

    fast = False
    if fast:
        header = "%mem=1000MB\n# AM1 SP"
        gaussian = (header + "\n\n" + smi + "\n" + 
                             "\n".join(mol.write("gau").split("\n")[3:]))
    else:
        gaussian = (header + "\n\n" + smi + "\n" + 
            "\n".join(mol.write("gau").replace("0  3\n", "0  1\n").split("\n")[3:])
            + header_b) % (idx, idx)

    output = open(os.path.join("gaussian", "%s.gjf" % idx), "w")
    output.write(gaussian)
    output.close()
