import os
import random
import pdb
import json
import pybel
import numpy as np

import Efficiency as effmod

# The following is taken from
#   http://wiki.python.org/moin/PythonDecoratorLibrary#Memoize


class memoized(object):
    """Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        try:
            return self.cache[args]
        except KeyError:
            value = self.func(*args)
            self.cache[args] = value
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)

    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__

    def __get__(self, obj, objtype):
        """Support instance methods."""
        return functools.partial(self.__call__, obj)

@memoized
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
    #newdimer = pybel.readstring("smi", reverse+reverse)
    #a, b = dimer.write("can"), newdimer.write("can")
    #assert a == b, "\n %s is not the same as\n %s" % (a, b)
    
    return reverse

@memoized
def issym(smiles):
    """Is this monomer symmetric?
    """
    rev = getreversed(smiles)
    return (pybel.readstring("smi", smiles+"Br").write("can") ==
            pybel.readstring("smi", rev + "Br").write("can"))

def _readdatafiles():
    relpath = os.sep.join(__file__.split(os.sep)[:-1])

    f = open(os.path.join(relpath, "directions_sequences.json"))
    # Need to remove the unicode from the json or else Open Babel won't like it
    combinations = dict([(float(x), [[str(t) for t in z] for z in y]) for x, y in json.load(f).iteritems()])
    f.close()
    return combinations

combinations = _readdatafiles()

def molname_to_repr(molname, length):
    """Convert a molecule 'name' to representation of the
    composition of the polymer

    >>> molname = "c(s1)c(OCO2)c2c1~C=CN=N_1_1"
    >>> myrepr = molname_to_repr(molname, 4)
    >>> print myrepr
    The base dimer is:
      C=CN=N attached to c(s1)c(OCO2)c2c1
    The polymer is formed from the base dimer as follows:
      + -
    The corresponding SMILES string is thus:
      C=CN=Nc(s1)c(OCO2)c2c1c(s1)c(OCO2)c2c1N=NC=C
    """
    parts, d1, d2 = molname.split("_")
    b = parts.split("~") # monos

    result = []
    comb = combinations[length][int(d1)][int(d2)]
    chomps = []
    for i in range(0, len(comb), 4):
        chomps.append(comb[i:i+4])
    base_dimer = ["", ""]
    for i in range(2):
        if chomps[i] == "(qu)":
            base_dimer[i] = b[0]
        elif chomps[i] == "(uq)":
            base_dimer[i] = "%s (reversed as %s) " % (b[0], getreversed(b[0]))
        elif chomps[i] == "(di)":
            base_dimer[i] = b[1]
        elif chomps[i] == "(id)":
            base_dimer[i] = "%s (reversed as %s) " % (b[1], getreversed(b[1]))
    result.append("The base dimer is:")
    result.append("  %s attached to %s" % (base_dimer[0], base_dimer[1]))

    startchomp = "".join(chomps[0] + chomps[1])
    myd = []
    for i in range(0, len(chomps), 2):
        if "".join(chomps[i] + chomps[i+1]) == startchomp:
            myd.append("+")
        else:
            myd.append("-")
            
    result.append("The polymer is formed from the base dimer as follows:")
    result.append("  " + " ".join(myd))
    result.append("The corresponding SMILES string is thus:\n  %s" % molname_to_smile(molname, length))
    return "\n".join(result)

def molname_to_smile(molname, length):
    """Convert a molecule 'name' to a molecule

    >>> molname = "c(s1)c(OCO2)c2c1~C=CN=N_1_0"
    >>> smile = molname_to_smile(molname, 4)
    >>> print smile
    C=CN=Nc(s1)c(OCO2)c2c1C=CN=Nc(s1)c(OCO2)c2c1
    """
    output = get_comb(molname, length)
    comb = output[0]
    b = output[1]
    print comb

    reversed = []
    for i in range(2):
        if issym(b[i]):
            reversed.append(b[i])
        else:
            reversed.append(getreversed(b[i]))
    smile = comb.replace("(qu)", b[0]).replace("(uq)", reversed[0]).replace(
                         "(di)", b[1]).replace("(id)", reversed[1])
    print smile
    return smile

def molname_to_mol(molname, length):
    smi = molname_to_smile(molname, length)
    return pybel.readstring("smi", smi)


def get_comb(molname, length):
    combinations = _readdatafiles()
    parts, d1, d2 = molname.split("_")
    comb = combinations[length][int(d1)][int(d2)]
    b = parts.split("~") # monos

    return comb, b


def polname(polymer, spacer="_"):
    b, d1, d2 = polymer
    ans = ["~".join(b), str(d1), str(d2)]
    return spacer.join(ans)

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

def dirindex(monos):
    idx = 0
    if issym(monos[0]):
        if issym(monos[1]):
            idx = 3
        else:
            idx = 1
    else:
        if issym(monos[1]):
            idx = 2

    if monos[0] == monos[1]:
        if issym(monos[0]):
            idx = 5
        else:
            idx = 4
    return idx

def alldirs(monos, length):
    """Return all possible directions"""
    idx = dirindex(monos)
    for comb in range(len(combinations[length][idx])):
        yield idx, comb

def randomdirs(monos, length):
    """Create random directions"""
    idx = dirindex(monos)
    comb = random.randint(0, len(combinations[length][idx]) - 1)
    return idx, comb

def besttrans(etens, etoscs, return_osc=False):
    """Return the scaling factor and energy of the best transition"""
    idx = 0
    found = -1
    max = 0
    # Finds largest etosc
    for idx in range(0, len(etens)):
        if etoscs[idx] >= 1.0 and found==-1:
            found = idx
        if etoscs[idx] > etoscs[max]:
            max = idx
    # Pick the etens that goes with the maximum oscillator strength
    scale = 1.0
    chosenidx = found
    if found == -1:
        # Use the strongest oscillator, and scale by its strength
        chosenidx = max
        scale = etoscs[chosenidx]

    if not return_osc:
        return scale, etens[chosenidx]
    else:
        return scale, etens[chosenidx], etoscs[chosenidx]

def besttrans_revised(etens, etoscs, return_osc=False):
    """Return the scaling factor and energy of the best transition"""
    # Pick the largest oscillator strength
    etoscs_max = max(etoscs)
    # Find the transition that corresponds to the largest oscillator strength
    etens_max_etoscs = etens[etoscs.index(max(etoscs))]

    # The scaling factor is the largest oscillator strength with a maximum value of 1
    scale = 1.0
    if etoscs_max <= 1.0:
        scale = max(etoscs)

    if not return_osc:
        return scale, etens_max_etoscs

    else:
        return scale, etens_max_etoscs, etoscs_max


def getHplusBG(json):
    """Get the HOMO and lowest energy significant transition
    given the JSON from the sqlite database"""
    homo, lumo, etens, etoscs = json.loads(json)
    scale, trans = besttrans(etens, etoscs)
    return homo, trans

class ScoreCalculator(object):
    def __init__(self):
        self.eff = effmod.Efficiency()
    def getscore(self, json, cutoff=True):
        homo, lumo, etens, etoscs = json.loads(json)
        scale, trans = besttrans(etens, etoscs)
        score = scale * self.eff.efficiency(homo, trans, -4.61, cutoff=cutoff)
        return score
    def getdistance(self, json):
        homo, lumo, etens, etoscs = json.loads(json)
        scale, trans = besttrans(etens, etoscs)
        penalty = 1.0 - scale
        distance = math.sqrt((homo-(-5.70))**2 + (trans-1.39)**2)
        score = distance + penalty
        return score

if __name__ == "__main__":

    import doctest
    doctest.testmod()
