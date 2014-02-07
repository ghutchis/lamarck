"""
maxdipole.py - Generate data for a Science paper

Usage: python maxdipole.py jobname [continue=N]

Two files are created, jobname.txt (the log file), and jobname.db,
an sqllite database.

To avoid overwriting previous results, if jobname.txt exists the
run will fail unless the continue=N option is provided (where N
is the step from which to continue).

jobname.db exists for the purpose of caching previous results. It
is only ever appended to.
"""
import os
import re
import sys
import pdb
import copy
import math
import time
import gzip
import random
import logging
import collections
import subprocess
import StringIO
from optparse import OptionParser

relpath = os.sep.join(__file__.split(os.sep)[:-1])

import json

from cclib.parser import ccopen, utils
import pybel

convert = 1.0 / utils.convertor(1, "eV", "cm-1")

from Efficiency import *
from Utils import *
from Admin import *

FAIL = -9999999

# These are the monomer used for the local search
chosen_monomers = {8: ['c(s1)c(OCCCO2)c2c1',
                       'c(s1)c(cccc2)c2c1',  'c(s1)cc(c12)[nH]c(c2s3)cc3',
                       'c(s1)cc(c12)Cc(c2s3)cc3', 'c(s1)cc(c12)c(=O)c(c2s3)cc3',
                       'c1oc(c2)c(c1)oc2', 'c(s1)cc(OC)c1',
                       'C(C1=C)=C(C=CC=C2)C2=C1', 'c(s1)c(ON=N2)c2c1',
                       'c(s1)c(OC=C2)c2c1', 'c(s1)c(NC=N2)c2c1',
                       'c(s1)c(SC=N2)c2c1', 'c(s1)c(OC=N2)c2c1',
                       'c(s1)c(nccn2)c2c1', 'c1oc(c2)c(c1)sc2',
                       'C(S1)=CC2=C1C=C(C2=O)', 'C1=CC(C2)=C(C1)C=C2',
                       'c(s1)c(cc(OC)c(OC)c2)c2c1', 'c(s1)c(OCCO2)c2c1',
                       'c(s1)c(NC=C2)c2c1', 'c(s1)c(ncnc2)c2c1'],
                   6: ['c(s1)c(ncnc2)c2c1',
                       'c(s1)c(NC=C2)c2c1',
                       'c(s1)c(nccn2)c2c1',
                       'c(s1)c(NC=N2)c2c1',
                       'c1oc(c2)c(c1)sc2',
                       'c1oc(c2)c(c1)oc2',
                       'c(s1)cc(c12)Cc(c2s3)cc3',
                       'c(s1)c(OC=C2)c2c1',
                       'C(S1)=CC2=C1C=C(C2=O)',
                       'c(s1)cc(c12)[nH]c(c2s3)cc3',
                       'c(s1)c(cc(OC)c(OC)c2)c2c1',
                       'c(s1)cc(c12)c(=O)c(c2s3)cc3',
                       'c(s1)cc(OC)c1',
                       'C1=CC(C2)=C(C1)C=C2',
                       'c(s1)c(OCCO2)c2c1',
                       'c(s1)c(OC=N2)c2c1',
                       'c(s1)c(OCCCO2)c2c1',
                       'c([nH]1)ccc1C=C',
                       'c(s1)c(cccc2)c2c1',
                       'c(s1)c(SC=N2)c2c1',
                       'c(s1)c(ON=N2)c2c1',
                       'c(s1)c(OCCS2)c2c1',
                       'c(s1)c(c(F)c(F)c(F)c(F)2)c2c1',
                       'c(s1)c(OC=CO2)c2c1',
                       'C(C1=C)=C(C=CC=C2)C2=C1',
                       'c1sc(c2)c(c1)sc2',
                       'c(s1)c(SN=N2)c2c1',
                       'c(s1)cc(O)c1',
                       'C(C1=O)=C(C=CC=C2)C2=C1',
                       'c1[nH]c(c2)c(c1)[nH]c2',
                       'c(s1)c(OCO2)c2c1',
                       'c(o1)ccc1',
                       'c(s1)c(cnnc2)c2c1',
                       'c([nH]1)ccc1C#C',
                       'c(s1)c(NCN2)c2c1',
                       'c(s1)ccc1C=C',
                       'c1sc(cc2c3)c(c1)cc2sc3',
                       'c(s1)c(cc(F)c(F)c2)c2c1',
                       'c(s1)cc(C=C)c1',
                       'c(s1)cc(c12)C(=C(C#N)C#N)c(c2s3)cc3',
                       'c([nH]1)cnc1',
                       'c(s1)cc(S)c1',
                       'c(s1)c(N(=O)=O)c(N)c1',
                       'c(s1)c(SCS2)c2c1',
                       'C(C1=O)=CC=C1',
                       'c(c(non1)c12)ccc2',
                       'C=CN=N',
                       'c(cc1)cc(c12)C(=O)c(c2c3)ccc3',
                       'c(s1)c(cc(C#N)c(C#N)c2)c2c1',
                       'c(s1)c(C(=O)OC2(=O))c2c1']
                   [:29]}


def exitwitherror():
    sys.exit(__doc__)


def createAllCombinations(mymonos, length):
    # Make all possible polymers of a particular length 
    dimerunits = []
    for i in range(len(mymonos)):
        for j in range(i, len(mymonos)):
            monos = [mymonos[i], mymonos[j]]
            if monos[1] > monos[0]: # In alphabetical order
                monos = [monos[1], monos[0]]
            ndir = 0
            for directions in alldirs(monos, length):
                dimerunits.append((monos, directions[0], directions[1]))
    return dimerunits


class GA(object):
    def __init__(self, admin, length, Nchromos, R, simmatrix, objectivefn,
                 logmessage=""):
        self.admin = admin
        self.N = Nchromos
        self.R = R # the number of nbrs
        self.simmatrix = simmatrix
        self.monomers = sorted(self.simmatrix.keys())
        self.length = length
        self.objectivefn = objectivefn
        self.gen = 0
        self.initscorefn()
        if logmessage:
            self.log(logmessage)

    def initscorefn(self):
        self.efficiency = Efficiency()
        self.efficiency.unittest()

    def log(self, msg):
        self.admin.log.write(msg + "\n")
        print msg

    def initpop(self):
        self.log("\tInitialising population")
        # Make self.N polymers of length self.length
        dimerunits = []
        for i in range(self.N):
            dimerunit = [[], "", ""]
            monos = [random.choice(self.monomers) for j in range(2)]
            if monos[1] > monos[0]: # In alphabetical order
                monos = [monos[1], monos[0]]
            directions = randomdirs(monos, self.length)
            dimerunits.append((monos, directions[0], directions[1]))
        self.pop = dimerunits
        self.logpop(self.pop)

    def initallpop(self, chosen_monos=None):
        if chosen_monos is None:
            chosen_monos = self.monomers
        self.log("\tInitialising population")
        self.N = 0

        dimerunits = createAllCombinations(chosen_monos, self.length)

        self.log("\tTotal size of potential pop is %d" % len(dimerunits))
        newunits = []
        for dimerunit in dimerunits:
            data = self.admin.getdata(polname(dimerunit))
            if not data:
                newunits.append(dimerunit)
        self.pop = newunits
        self.log("\tTotal size of uncalc pop is %d" % len(newunits))
        self.logpop(self.pop)

    def logpop(self, pop):
        self.log("Population of size %d" % len(pop))
        for i, pol in enumerate(pop):
##            self.log("%d: %s %s %s" % (i, pol[0], pol[1], pol[2]))
            self.log(polname(pol))

    def loggjf(self):
        self.log("%d GJF files created" % len(self.gjfs))
        for j, x in enumerate(self.gjforder):
            self.log("GJF %d: polymer numbers %s" % (j, self.gjfs[x]))

    def logfitness(self, pop, text):
        self.log("\t%s population fitness" % text)
        for j, x in enumerate(sorted(pop, key = lambda x: self.getscore(polname(x)), reverse=True)):
            score, logtext = self.getscore(polname(x), log=True)
            if score is not None:
                self.log("%d: %s with %.3f %s" % (j, polname(x), score, logtext))
            else:
                self.log("%d: %s with FAIL" % (j, polname(x)))

    def makeGJF(self, pop, length):
        self.log("\tCreating txt files")
        if not os.path.isdir("gaussian"):
            os.mkdir("gaussian")
        self.gjfs = {}

        for i, x in enumerate(pop):
            if self.getscore(polname(x)) == None:
                self.gjfs.setdefault(polname(x), []).append(i)

                print "showing the pol: %s" % polname(x)

        # Sort the gjfs by molecular weight
        self.gjforder = sorted(self.gjfs.keys(), key=lambda x: molname_to_mol(str(x), self.length).molwt, reverse=True)

        for idx, smi in enumerate(self.gjforder):
            mol = molname_to_mol(smi, length)
            mol.make3D()
            globalopt(mol)

            header = "%%nproc=1\n%%mem=1GB\n%%Chk=%s.chk\n#T PM6 OPT"
            header_b = """
--Link1--
%%nproc=1
%%mem=1GB
%%Chk=%s.chk
%%NoSave
# Geom=AllCheck ZINDO(NStates=15,Singlets)
"""
            gaussian = (header + "\n\n" + smi + "\n"
                + "\n".join(mol.write("gau").replace("0  3\n", "0  1\n").split("\n")[3:])
                + header_b) % (idx, idx)
            output = open(os.path.join("gaussian", "%s.gjf" % idx), "w")
            output.write(gaussian)
            output.close()
            print "finished creating %s.gjf" % idx

    def runGaussian(self):
        CPUS_PER_NODE = 4

        if len(self.gjfs) > 0:
            if os.path.isfile(os.path.join("gaussian", "end.txt")):
                self.log("\nFound end.txt. Finishing")
                sys.exit(0)
            self.log("\tRunning Gaussian")
            output = open("tasks", "w")
            for i in range(len(self.gjfs)):
                output.write(
                    "echo -n \"%d \"; date; cd $PBS_O_WORKDIR; python %s/smi23D.py %d %d; cd gaussian; g09 < %d.gjf > %d.out; gzip %d.out; echo -n \"%d \"; date; rm -f %d.chk\n" % (i, relpath, i, self.length, i, i, i, i, i))

            output.close()
            template = open("template.sh", "r").read()
            time_per_job = 3 # 16 for 8mers, 11 for 6, 6 for 4, 3 for dimers

            if self.length >= 8 and len(self.gjfs) <= 64:
                # Schedule the long jobs first, and give those slow jobs
                # extra time
                scheme = "largevariation"
            else:
                scheme = "normal"

            if scheme == "largevariation":
                N_parallel_jobs = len(self.gjfs) * 5 / CPUS_PER_NODE
                if N_parallel_jobs == 0:
                    N_parallel_jobs = 1
                mins = time_per_job * 2
                hours = 0
                N_nodes = 1 + (N_parallel_jobs - 1) / CPUS_PER_NODE
            else: # Normal
                MAX = CPUS_PER_NODE * 8
                if len(self.gjfs) <= MAX:
                    M = 1
                    N_nodes = 1 + (len(self.gjfs) - 1) / CPUS_PER_NODE
                else:
                    M = 1 + (len(self.gjfs) - 1) / MAX
                    N_nodes = MAX / CPUS_PER_NODE
                walltime = int(M * time_per_job)
                mins = walltime % 60
                hours = walltime / 60
            # Workaround for Stokes (3 nodes not allowed)
            if N_nodes == 3:
                N_nodes= 4
            template = template.replace("REPLACENODES", str(N_nodes))
            template = template.replace("REPLACEHOUR", str(hours))
            template = template.replace("REPLACEMIN", str(mins))
            output = open("runwith1.sh", "w")
            print >> output, template
            output.close()
            for i in range(len(self.gjfs)):
                filename = os.path.join("gaussian", "%d.out.gz" % i)
                if os.path.isfile(filename):
                    os.remove(filename)
            qsub = subprocess.Popen(["qsub", "runwith1.sh"], stdout=subprocess.PIPE)
            stdout = qsub.stdout.read()
            self.log(stdout)
            pid = stdout.split(".")[0]
            stderr = ""
            while not stderr.strip():
                time.sleep(10)
                qstat = subprocess.Popen(["qstat", pid], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout = qstat.stdout.read()
                stderr = qstat.stderr.read()

    def extractcalcdata(self):
        if len(self.gjfs) > 0:
            self.log("\tExtracting data from log files")
            tostore = []
            for j, pname in enumerate(self.gjforder):
                mylogfile = os.path.join("gaussian", "%d.out.gz" % j)
                if not os.path.isfile(mylogfile):
                    continue
                text = gzip.open(mylogfile, "r").read()
                if text.find("Excitation energies and oscillator strength") < 0:
                    continue
                lines = iter(text.split("\n"))
                for line in lines:
                    if line.startswith(" #T PM6 OPT"):
                        line = lines.next()
                        line = lines.next()
                        line = lines.next()
                        break
                for line in lines:
                    if line.startswith(" Initial command"): break
                zindofile = list(lines)
                if len(zindofile) == 0:
                    # All the PM6 data is missing
                    continue
                with open("tmp.out", "w") as f:
                     f.write("\n".join(zindofile))
                logfile = ccopen("tmp.out")
                logfile.logger.setLevel(logging.ERROR)
                try:
                    data = logfile.parse()
                except AssertionError:
                    continue
                try:
                    lumo = data.moenergies[0][data.homos[0] + 1]
                    homo = data.moenergies[0][data.homos[0]]
                    etens = [x*convert for x in data.etenergies] # cm-1 to eV
                    etoscs = data.etoscs
                except:
                    continue
                if max(etens) <= 0:
                    continue

##                myjson = json.dumps([homo, lumo, etens, etoscs])
                myjson = json.dumps([homo, lumo, etens, etoscs,
                                    data.moenergies[0], data.homos[0]])
                tostore.append((pname, myjson))

            for pname, myjson in tostore:
                self.admin.storedata(pname, self.gen, myjson)

    def getscore(self, polname, log=False):
        data = self.admin.getdata(polname)
        logtext = ""
        if not data:
            if log:
                return None, None
            else:
                return None
        gen, myjson = data
        jsondata = json.loads(myjson)
        if len(jsondata) == 4:
            homo, lumo, etens, etoscs = jsondata
        else:
            homo, lumo, etens, etoscs, moenergies, homo_idx = jsondata
        scale, trans = besttrans(etens, etoscs)
        if scale < 1.0:
            logtext += "Os=%.1f" % (scale*100,)

        if self.objectivefn == "eff":
            score = scale * self.efficiency.efficiency(homo, trans, -4.61)
        elif self.objectivefn == "distance":
            penalty = 1.0 - scale
            distance = math.sqrt((homo-(-5.70))**2 + (trans-1.39)**2)
            score = distance + penalty
            score = -score # We are finding the maximum
        if log:
            return score, logtext
        else:
            return score

    def makechildren(self, moverandomly=False):
        """
It should be possible for a single monomer to mutate

The mutations should always allow the exploration of local space
"""
        self.gen += 1
        scores = []
        for chromo in self.pop:
            x = polname(chromo)
            scores.append((chromo, self.getscore(x)))

        poolsize = self.N / 5
        pool = []
        for i in range(poolsize):
            tournament = random.sample(scores, 3)
            tournament.sort(reverse=True, key=lambda x:x[1])
            select = tournament[0]
            scores.remove(select)
            pool.append(select[0])

        self.children = []
        while len(self.children) < self.N:
            # Crossover to make two children
            x = copy.deepcopy(random.choice(pool))
            y = copy.deepcopy(random.choice(pool))
            children = [[x[0][0], y[0][1]], [x[0][1], y[0][0]]]
            for child in children:
                newchild = [child[0], child[1]]
                # Mutate backbone
                for i, mon in enumerate(child):
                    if random.random() > 0.25:
                        if moverandomly:
                            newchild[i] = random.choice(self.monomers)
                        else:
                            newchild[i] = random.choice(
                                           self.simmatrix[mon][:self.R])

                if newchild[1] > newchild[0]: # Alphabetical order
                    newchild = [newchild[1], newchild[0]]

                # Create random dirs
                directions = randomdirs(newchild, self.length)
                fullchild = (newchild, directions[0], directions[1])
                # Don't add a duplicate
                if fullchild not in self.pop + self.children:
                    self.children.append(fullchild)

    def nextgen(self):
        self.pop.sort(key=lambda x: self.getscore(polname(x)),
                      reverse=True)
        self.children.sort(key=lambda x: self.getscore(polname(x)),
                           reverse=True)

        self.pop = self.pop[:self.N/2] + self.children[:self.N/2]


def testdb(admin):
    pdb.set_trace()
    sys.exit(0)


def doGA(admin, length, Nchromos, nbrs, simmatrix,
         moverandomly=False, no_gaussian=False,
         objectivefn="eff", logmessage=""):
    """If no_gaussian is True, this will cause the algorithm not
    to generate and run any new Gaussian jobs. Any missing values in the
    database will be regarded as FAILs"""
    ga = GA(admin, length, Nchromos, nbrs, simmatrix, objectivefn, logmessage)
    ga.initpop()
    if not no_gaussian:
        ga.makeGJF(ga.pop, length)
        ga.runGaussian()
        ga.extractcalcdata()
    ga.logfitness(ga.pop, "Initial")
    for i in range(100):
        ga.makechildren(moverandomly)
        if not no_gaussian:
            ga.makeGJF(ga.children, length)
            ga.runGaussian()
            ga.extractcalcdata()
        ga.logfitness(ga.children, "Children")
        ga.nextgen()
        ga.logfitness(ga.pop, "Next generation")


def doExhaustive(admin, length, simmatrix, initallpop=False):
    ga = GA(admin, length, None, None, simmatrix, None)
    use_just_chosen_monomers = False
    if use_just_chosen_monomers:
        monos_to_use = chosen_monomers[length]
    else:
        monos_to_use = sorted(simmatrix.keys())

    print "The number of monomers is %d" % len(monos_to_use)
    if initallpop:
        myranges = [1, 5, 29, len(monos_to_use)]
        for x in myranges:
            print "Working on up to %d monomers" % x
            ga.initallpop(monos_to_use[:x])
            ga.makeGJF(ga.pop, length)
            ga.runGaussian()
            ga.extractcalcdata()
    else:
        for i in range(0, len(monos_to_use)):
            print "== Working on %d, %s ==" % (i, monos_to_use[i])
            ga.initpop(i)
            ga.makeGJF(ga.pop, length)
            ga.runGaussian()
            ga.extractcalcdata()


def test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    parser = OptionParser()
    parser.set_defaults(dotest=False, moverandomly=False, no_gaussian=False,
                        seed=1, database=None, Nchromos=16, length=2,
                        nbrs=None, simmatrix="", objectivefn="eff",
                        exhaustive=False, initall=False)
    parser.add_option("--initall", dest="initall",
                      help="in an exhaustive search, initialise all of the polymers",
                      action="store_true")
    parser.add_option("-x", "--exhaustive", dest="exhaustive",
                      help="do an exhaustive search",
                      action="store_true")
    parser.add_option("-s", "--seed", dest="seed", help="set the random seed",
                      type="int", metavar="SEED")
    parser.add_option("-d", "--database", dest="database",
                      help="set the database name",
                      metavar="DATABASE")
    parser.add_option("-t", "--test", dest="dotest", action="store_true",
                      help="run the test")
    parser.add_option("--no-gaussian", dest="no_gaussian",
                      action="store_true", help="don't run any Gaussian jobs")
    parser.add_option("-r", "--random", dest="moverandomly",
                      action="store_true",
                      help="chose the next generation randomly")
    parser.add_option("-l", "--length", dest="length", type="int",
                      help="set the polymer length", metavar="LENGTH")
    parser.add_option("-N", dest="Nchromos", type="int",
                      help="set the number of chromosomes", metavar="SIZE")
    parser.add_option("-R", "--nbrs", dest="nbrs", type="int", metavar="NBRS",
                      help="set the number of nearest neighbours")
    parser.add_option("--matrix", dest="simmatrix",
                      metavar="FILENAME",
                      help="which similarity matrix to use")
    parser.add_option("--objective", dest="objectivefn", type="choice",
                      choices=["eff", "distance"], metavar="OBJ_FN",
                      help="which objective function (eff or distance)")

    (options, args) = parser.parse_args()
    if not options.exhaustive:
        if not options.database:
            parser.error("no database specified")
        if not options.nbrs:
            parser.error("you need to specify the number of neighbours")

    if not options.simmatrix or not os.path.isfile(options.simmatrix):
        parser.error("you need to specify the similarity matrix")
    else:
        with open(options.simmatrix, "r") as f:
            simmatrix = json.load(f)
            # Remove unicode from JSON
            tmp = {}
            for x, y in simmatrix.iteritems():
                tmp[str(x)] = map(str, y)
            simmatrix = tmp

    admin = Admin(options.database)
    if options.dotest:
        testdb(admin)
    else:
        if not options.exhaustive:
            random.seed(options.seed)
            doGA(admin, options.length, options.Nchromos, options.nbrs,
                 simmatrix, logmessage = "Chosen options: %s" % str(options),
                 moverandomly=options.moverandomly,
                 no_gaussian=options.no_gaussian,
                 objectivefn=options.objectivefn)
        else:
            doExhaustive(admin, options.length,
                         simmatrix, initallpop=options.initall)
