"""
"""
import os
import sys
import pdb

import pybel
import simplejson

from Utils import molname_to_mol, globalopt

def exitwitherror():
    sys.exit(__doc__)

header = "%%nproc=4\n%%mem=512MB\n%%Chk=%s.chk\n#T PM6 OPT"
header_b = """
--Link1--
%%nprocshared=4
%%Chk=%s.chk
%%NoSave
#t Geom=AllCheck B3LYP/6-31G* SP Prop=FitCharge
-"""

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

