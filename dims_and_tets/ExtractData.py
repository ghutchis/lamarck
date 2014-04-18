import os
import sys
import gzip
import glob
import logging
import StringIO

from cclib.parser import ccopen, utils
import pybel

def extractdata(folder):
    smiles = [x.rstrip() for x in open(os.path.join(folder, os.path.basename(folder) + ".txt"), "r").readlines()]
    print smiles
    archivefile = open(os.path.join(folder, "zindo.txt"), "w")
    print >> archivefile, "\t".join(["File ID", "SMILES", "HOMO (eV)", "LUMO (eV)",
               "Trans (eV)", "Osc", "..."])

    getnum = lambda x: int(x.split("/")[1].split(".")[0])
    homos = []
    lumos = []
    trans = []
    convert = 1.0 / utils.convertor(1, "eV", "cm-1")
    for filename in sorted(glob.glob(os.path.join(folder, "*.gz")),
                           key=getnum):

        number = getnum(filename)
        smile = smiles[number]
        text = gzip.open(filename, "r").read()

        if text.find("Excitation energies and oscillator strength") < 0:
            continue
        lines = iter(text.split("\n"))

        for line in lines:
            if line.startswith(" Initial command"): break
        text = StringIO.StringIO("\n".join(list(lines)))

        logfile = ccopen(text)
        logfile.logger.setLevel(logging.ERROR)
        data = logfile.parse()

        #assert(len(data.homos) == 1)
        smiles.append(smile)
        homo = data.homos[0]
        homos.append(data.moenergies[0][homo])
        lumos.append(data.moenergies[0][homo + 1])
        trans.append(zip(data.etenergies, data.etoscs))

        archivefile.write("%d\t%s\t%f\t%f" % (number, smile, homos[-1], lumos[-1]))
        for x in trans[-1]:
            archivefile.write("\t%f\t%f" % (x[0] * convert, x[1]))
        archivefile.write("\n")
##        print >> open("tmp.txt", "w"), text.getvalue()
        if smile != "c(s1)c(SN=N2)c2c1c(s1)c(SN=N2)c2c1c(s1)c(SN=N2)c2c1c(s1)c(SN=N2)c2c1c(s1)c(SN=N2)c2c1c(s1)c(SN=N2)c2c1":
            mol = pybel.readstring('g09', text.getvalue())
            mol.write("xyz", os.path.join(folder, "%d.xyz" % number), overwrite=True)

    print "%s: Created zindo.txt, plus various xyz files." % folder

def showhelp():
    print """
  Usage: extractdata.py folder

         where folder contains the results of series of
         calculation on homopolymers."""
    sys.exit(1)

if __name__ == "__main__":
##    for folder in ['monolen2', 'monolen4', 'monolen6', 'monolen8']:
##        extractdata(folder)
    if len(sys.argv) != 2:
        showhelp()
    folder = sys.argv[1]
    if not os.path.isdir(folder):
        showhelp()
    extractdata(folder)