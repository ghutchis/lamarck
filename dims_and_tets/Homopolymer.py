import os
import sys
import pybel
from optparse import OptionParser

sys.path.append(os.path.join("..", "geneticAl"))
#from Utils import globalopt

header = "%%nproc=1\n%%mem=1GB\n%%Chk=%d.chk\n#T PM6 OPT=(MaxCycles=500)"
header_b = """
--Link1--
%%nproc=1
%%mem=1GB
%%Chk=%d.chk
%%NoSave
# Geom=AllCheck ZINDO(NStates=15,Singlets)
"""
def globalopt(mol, debug=False, fast=False):
    pybel._builder.Build(mol.OBMol)
    mol.addh()
    if debug:
        ff.GetCoordinates(mol.OBMol)
        mol.write("sdf", "1.sdf", overwrite=True)

    ff = pybel._forcefields["mmff94"]
    success = ff.Setup(mol.OBMol)
    if not success:
        ff = pybel._forcefields["uff"]
        success = ff.Setup(mol.OBMol)
        if not success:
            sys.exit("Cannot set up forcefield")

    if fast:
        ff.SteepestDescent(50, 1.0e-3)
    else:
        ff.SteepestDescent(500, 1.0e-4)
    if debug:
        ff.GetCoordinates(mol.OBMol)
        mol.write("sdf", "2.sdf", overwrite=True)

    if fast:
        ff.WeightedRotorSearch(20, 5)
    else:
        ff.WeightedRotorSearch(100, 20)
    if debug:
        ff.GetCoordinates(mol.OBMol)
        mol.write("sdf", "3.sdf", overwrite=True)

    if fast:
        ff.ConjugateGradients(50, 1.0e-4)
    else:
        ff.ConjugateGradients(500, 1.0e-6)
    ff.GetCoordinates(mol.OBMol)
    if debug:
        mol.write("sdf", "4.sdf", overwrite=True)

def getmonomers(filename, debug=False):
    """
    >>> monos = getmonomers()
    The number of monomers is 142
    """
    inputfile = filename
    exclude = ["se", "p", "b", "C1=CC(C2=O)=C(C1=O)C=C2",
               "C(=C1)C=C1"]
    excluded = []
    monos = []
    for smile in open(inputfile, "r"):
        if any(smile.find(x) >= 0 for x in exclude):
            excluded.append(smile)
        else:
            monos.append(smile.rstrip())
    oldlen = len(monos)
    monos = set(monos)
    print excluded
    print monos

    if debug:
        print "".join(excluded)
        print "%d duplicates removed" % (oldlen - len(monos))
    
    print "The number of monomers is", len(monos)
    return list(monos)
        
def makefiles(filename, folder, length, repunit = 1):
    monos = getmonomers(filename)

    smiles = []
    for i, smile in enumerate(monos):
        if repunit == 1:
            smiles.append( smile*length )
        elif repunit == 2:
            smiles.extend( ["%s%s" % (smile, x) * (length / repunit)
                            for x in monos[i+1:]])

    info = open(os.path.join(folder, os.path.basename(folder) + ".txt"), "w")
    sdf = pybel.Outputfile("sdf", os.path.join(folder, os.path.basename(folder) + ".sdf"), overwrite=True)
    for i, smile in enumerate(smiles):
        print i, str(smile)
        print >> info, i, str(smile)
        mol = pybel.readstring("smi", smile)
        globalopt(mol)
        gaussian = (header + "\n\n" + smile + "\n" + 
                    "\n".join(mol.write("gau").replace("0  3\n", "0  1\n").split("\n")[3:])
                    + header_b) % (i, i)
        with open(os.path.join(folder, "%d.com" % i), "w") as output:
            output.write(gaussian)
        mol.title = str(i)
        sdf.write(mol)
    info.close()
    sdf.close()

def test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":

    parser = OptionParser()
    parser.set_defaults(monomerfilename=None)
    parser.add_option("-f", "--monomerfile", dest="monomerfilename",
                      help="file containing the monomers",
                      metavar = "FILENAME")
    parser.add_option("-d", "--directory", dest="folder",
                      help="existing directory to be used for gaussian input files",
                      metavar = "DIRECTORY")
    parser.add_option("-l", "--length", dest="length", type="int",
                      help="set the polymer length", metavar="LENGTH")

    (options, args) = parser.parse_args()
    if not options.monomerfilename:
        parser.error("no monomer filename specified")
    if not options.folder or not os.path.isdir(options.folder):
        parser.error("no existing directory specified")

    print options
    makefiles(options.monomerfilename, options.folder, options.length)

