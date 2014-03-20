import os
import sys
import pybel
from optparse import OptionParser

sys.path.append(os.path.join("..", "geneticAl"))
from Utils import globalopt

header = "%%nproc=1\n%%mem=1GB\n%%Chk=%d.chk\n#T PM6 OPT"
header_b = """
--Link1--
%%nproc=1
%%mem=1GB
%%Chk=%d.chk
%%NoSave
# Geom=AllCheck ZINDO(NStates=15,Singlets)
"""

def getmonomers(filename, debug=False):
    """
    >>> monos = getmonomers()
    The number of monomers is 142
    """
    inputfile = filename
    exclude = ["se", "p", "Si", "b", "C1=CC(C2=O)=C(C1=O)C=C2",
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

    info = open(os.path.join(folder, folder + ".txt"), "w")
    sdf = pybel.Outputfile("sdf", os.path.join(folder, folder + ".sdf"), overwrite=True)
    for i, smile in enumerate(smiles):
        print i, str(smile)
        print >> info, str(smile)
        mol = pybel.readstring("smi", smile)
        globalopt(mol)
        gaussian = (header + "\n\n" + smile + "\n" +
                    "\n".join(mol.write("gau").replace("0  3\n", "0  1\n").split("\n")[3:])
                    + header_b) % (i, i)
        with open(os.path.join(folder, "%d.gjf" % i), "w") as output:
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