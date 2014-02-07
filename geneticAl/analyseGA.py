import os
import sys
import pdb
import math
import pylab
import matplotlib as mpl
import numpy as np
import pickle
import collections
try:
    import simplejson
except ImportError:
    import json as simplejson

sys.path.append(os.path.join("..", "dims_and_tets"))
import solar
    
from Efficiency import *
from Utils import *
from Admin import *

import geneticAl

chosentetramers = set(x.split()[0] for x in open(os.path.join(
                      "..", "dims_and_tets", "best_alltetramers.txt"), "r"))
assert len(chosentetramers) == 101
besttetramers = set(x.split()[0] for x in open(os.path.join(
                      "..", "dims_and_tets", "best_alltetramers.txt"), "r")
                    if float(x.split()[-1]) > 9.23)
assert len(besttetramers) == 10 # Top 10

class GA_Results(object):
    def __init__(self, filename, database, length):
        self.filename = filename
        self.database = database
        self.length = length
        self.eff = Efficiency()
        self.admin = Admin(self.database)
        self.readlogfile()

    def readlogfile(self):
        log = open("%s" % self.filename, "r")
        line = log.next()
        while not line.startswith("Population of size"):
            line = log.next()
        self.size = int(line.rstrip().split()[-1])
            
        for line in log:
            if line.find("Initial population fitness") >= 0:
                break
        gen = 0
        self.vals = []
        self.xy = []
        vals = [[], [], [], []]
        pol = set()
        allpop = []
        for gen in range(100):
            plotxy = []
            data = []
            for i in range(self.size):
                line = log.next()
                broken = line.rstrip().split()
                if gen == 0:
                    pol.add(molname_to_smile(broken[1], self.length))
                if not broken[3] == "FAIL":
                    gen, json = self.admin.getdata(broken[1])                    
                    lookup = getHplusBG(json)
                    data.append( (broken[1], float(broken[3]),
                                  lookup[0], lookup[1]) )
            allpop.append( [(x[0], x[1]) for x in data] )

            vals[0].append(np.mean([x[1] for x in data]))
            vals[1].append(np.mean([x[1] for x in data[:(len(data)/2)]]))
            vals[2].append(data[0][1])
            vals[3].append(len(pol & chosentetramers))
                        
            plotxy = ([x[2] for x in data], [x[3] for x in data])
            self.xy.append(plotxy)

            for line in log:
                if line.find("Children population fitness") >= 0:
                    break
            # Either at end of file or start of Children
            if line.find("Children population fitness") >= 0:
                for i in range(self.size):
                    broken = log.next().rstrip().split()
                    pol.add(broken[1])
                
                for line in log:
                    if line.find("Next generation population fitness") >= 0:
                        break
        self.vals = vals
        self.pols = pol
        self.allpop = allpop
        
    def histo(self):
        data = []
        for x in self.admin.getalldata():
            e = self.getscore(x[2])
            if e>0:
                data.append(e)
        pylab.hist(data, bins=range(0, 14))
        pylab.show()

    def plot(self, prefix="noel"):
        names = ['mean', 'mean16', 'max']
        for i in range(3):
            pylab.plot(range(101), self.vals[i])
        pylab.savefig("pictures/prog%s.png" % prefix)
        pylab.clf()

    def progress(self, prefix="noel"):
        for N in range(20):
            solar.plot()            
            pylab.plot(self.xy[N][1], self.xy[N][0], ".")
            pylab.xlim(0.5, 4.5)
            pylab.ylim(-8.5, -4.5)
            pylab.savefig("pictures/%s%02d.png" % (prefix, N))
            pylab.clf()
##        pylab.show()

def novel_polymers():
    a = []
    b = []
    for i in [8, 16, 32, 64, 128, 256]:
        res = GA_Results("dimer_ga/%d_1.txt" % i, "alldimers")
        print i
        a.append(len(res.pols))
        b.append((len(res.pols)*100) / float(i*101)) # %
    print a
    print b

def investigate_params():
    dbs = {2:'alldimers', 4:'alltetramers'}
    base = os.path.join("paramsweep", "output")
    for length in [4]:
        pols = dbs[length]
        db  = os.path.join("..", "dims_and_tets", pols)
        for i in range(10):
            res = GA_Results(os.path.join(base,
                    "%s_CHR64_NBR7_MAT4_distance_SEED%d.txt" % (pols, i)),
                             db, length)
            pylab.plot(range(100), res.vals[3], label="%d" % i)
            print len(res.pols), len(res.pols & chosentetramers), len(res.pols & besttetramers)
        pylab.legend(loc="lower right")
        pylab.xlim(0, 100)
        pylab.savefig(os.path.join("pictures", "%s_64_repeats.png" % pols))
        pylab.clf()
        
##        for nbr in [4, 7, 11, 15, 19, 70]:
##            res = GA_Results("dimer_ga/64_1_NBR%d_SIM4.txt" % nbr, "alldimers")
##            pylab.plot(range(100), res.vals[2], label="%d" % nbr)
##        pylab.legend(loc="lower right")
##    ##    pylab.xlim(0, 40)
##        pylab.savefig("pictures/effectNnbrs_64_1_SIM4.png")
##
##        for matrix in [2, 4]:
##            pylab.clf()
##            for i in [8, 16, 32, 64, 128, 256]:
##                res = GA_Results("dimer_ga/%d_1_NBR7_SIM%d.txt" % (i, matrix), "alldimers")
##                pylab.plot(range(100), res.vals[2], label="%d" % i)
##            pylab.legend(loc="lower right")
##            pylab.ylim(0, 10)
##            pylab.savefig("pictures/effectNchromos_SIM%d.png" % matrix)    

def within(distance, (x1, y1), (x2, y2)):
    d = math.sqrt((x1-x2)**2 +(y1-y2)**2)
    return d < distance

def closetoline(distance, (x0, y0), (x1, y1), (x2, y2)):
    d = abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)) / math.sqrt(
                (x2-x1)**2 + (y2-y1)**2)
    return d < distance

def plotGA(filename, db, zoomlevel=1, save=False):
    pcbm_lumo = -4.61
    solar.plot(pcbm_lumo)    
    res = GA_Results(filename, db)
    admin = Admin(db)
    data = admin.getalldata()
    pdata = [getHplusBG(x[2]) for x in data]
    pylab.plot([x[1] for x in pdata], [x[0] for x in pdata], ".")

    for i in range(100):
        x, y = res.xy[i][1], res.xy[i][0]
        pylab.plot(x, y, "r.")

    limits = {2: ((1.15, 2.2), (-7.1, -5.4)),
                      1: ((0.5, 3.2), (-8.1, -4.8)),
                      0: ((0.5, 5.3), (-11.5, -3.7)),
                      }
    zoomlevel = 2
    pylab.xlim(*limits[zoomlevel][0])
    pylab.ylim(*limits[zoomlevel][1])
    if not save:
        pylab.show()
    else:
        base = os.path.basename(filename).split(".")[0]
        pylab.savefig(os.path.join("pictures", "%s.png" % base))
        pylab.clf()

def compare_objective_fns():
    """Decide which objective function is best

    Answer: distance is better than eff
    """
    for pols in ["alltetramers", "alldimers"]:
        for seed in range(10):
            plotGA(os.path.join("paramsweep", "output",
                                "%s_CHR64_NBR7_MAT4_distance_SEED%d.txt" % (pols, seed)),
                   os.path.join("..", "dims_and_tets", pols),
                   save=True)
            plotGA(os.path.join("paramsweep", "output",
                                "%s_CHR64_NBR7_MAT4_eff_SEED%d.txt" % (pols, seed)),
                   os.path.join("..", "dims_and_tets", pols),
                   save=True)

def compare_chromos():
    """Decide which number of chromos is best

    *** 8 ***
    [625.39999999999998, 9.0999999999999996, 0.90000000000000002]
    *** 16 ***
    [720.20000000000005, 30.399999999999999, 2.2999999999999998]
    *** 32 ***
    [1151.4000000000001, 40.200000000000003, 4.2000000000000002]
    *** 64 ***
    [2357.4000000000001, 58.700000000000003, 7.2000000000000002]
    *** 128 ***
    [3845.0, 74.900000000000006, 7.7999999999999998]
    *** 256 ***
    [5753.1999999999998, 81.0, 8.1999999999999993]
    *** 512 ***
    [10138.5, 92.0, 9.8000000000000007]
    """
    recalculate = False
    chromos = [8, 16, 32, 64, 128, 256, 512]
    if recalculate:
        alldata = []
        for chromo in chromos:
            print "***", chromo, "***"
            data = [[], [], []]
            for seed in range(10):
                filename = os.path.join("paramsweep", "output",
                                    "alltetramers_CHR%d_NBR7_MAT4_distance_SEED%d.txt" % (chromo, seed))
                db = os.path.join("..", "dims_and_tets", "alltetramers")
                res = GA_Results(filename, db, 4)
                result = len(res.pols), len(res.pols & chosentetramers), len(res.pols & besttetramers)
                for i in range(3):
                    data[i].append(result[i])
            print [pylab.mean(x) for x in data]
            alldata.append(data)
        with open("compare_chromos.pickle", "w") as f:
            pickle.dump(alldata, f)
    
    with open("compare_chromos.pickle", "r") as f:
        alldata = pickle.load(f)

    for i, x in enumerate(['Number of polymers', 'Number of chosen tetramers',
              'Number of top 10 most efficient tetramers']):
        pylab.boxplot([y[i] for y in alldata])
        pylab.xlabel("Number of chromosomes")
        pylab.ylabel(x)
        pylab.gca().set_xticklabels(chromos)
        pylab.savefig(os.path.join("pictures", "Nchromos_%d.png" % i))
        pylab.clf()
    pylab.plot(chromos, [pylab.mean(x[0]) for x in alldata])
    pylab.show()

def compare_nbrs():
    """Decide which number of nbrs is best

    *** 3 ***
    [1428.5, 49.100000000000001, 6.2999999999999998]
    *** 7 ***
    [2240.9000000000001, 56.799999999999997, 6.4000000000000004]
    *** 11 ***
    [3262.0, 66.900000000000006, 7.7999999999999998]
    *** 15 ***
    [4020.5, 66.299999999999997, 7.9000000000000004]
    *** 31 ***
    [5613.6999999999998, 70.700000000000003, 8.0]
    *** 66 ***
    [6858.1999999999998, 55.200000000000003, 6.5]
    *** 3 ***
    [1489.7, 40.700000000000003, 5.0999999999999996]
    *** 7 ***
    [2357.4000000000001, 58.700000000000003, 7.2000000000000002]
    *** 11 ***
    [3079.3000000000002, 66.700000000000003, 7.2000000000000002]
    *** 15 ***
    [3791.8000000000002, 62.799999999999997, 7.0999999999999996]
    *** 31 ***
    [5291.3999999999996, 68.799999999999997, 8.0999999999999996]
    *** 66 ***
    [6714.8999999999996, 60.600000000000001, 7.0999999999999996]
    """
    for matrix in [2, 4]:
        nbrs = [3, 7, 11, 15, 31, 66]

        alldata = []
        for nbr in nbrs:
            print "***", nbr, "***"
            data = [[], [], []]
            for seed in range(10):
                filename = os.path.join("paramsweep", "output",
                                        "alltetramers_CHR64_NBR%d_MAT%d_distance_SEED%d.txt" % (nbr, matrix, seed))
                db = os.path.join("..", "dims_and_tets", "alltetramers")
                res = GA_Results(filename, db, 4) # Check
                result = len(res.pols), len(res.pols & chosentetramers), len(res.pols & besttetramers)
                
                for i in range(3):
                    data[i].append(result[i])
            print [pylab.mean(x) for x in data]
            alldata.append(data)

        for i, x in enumerate(['Number of polymers', 'Number of chosen tetramers',
                  'Number of top 10 most efficient tetramers']):
            pylab.boxplot([y[i] for y in alldata])
            pylab.xlabel("Number of neighbours")
            pylab.ylabel(x)
            pylab.gca().set_xticklabels(nbrs)
            pylab.savefig(os.path.join("pictures", "Nnbrs_matrix%d_%d.png" % (matrix, i)))
            pylab.clf()

def find_best_monos():
    from Figure1 import best_monomers

    bestmonos, toponepc = best_monomers()
    print
    
    for seed in range(10):
        filename = os.path.join("paramsweep", "output",
                                "alltetramers_CHR64_NBR7_MAT4_distance_SEED%d.txt" % seed)
        db = os.path.join("..", "dims_and_tets", "alltetramers")
        res = GA_Results(filename, db, 4)
        pop = res.allpop

        data = {}

        for p in pop:
            data.update(dict(p))
        counts = collections.defaultdict(int)
        for x in data.keys():
            monos = x.split("_")[0].split("~")
            counts[monos[0]] += 1
            counts[monos[1]] += 1
##        pylab.hist(counts.values(), bins=range(0, 150, 2))
##        pylab.show()

        # If we generate all those dimers that involve the top 10 monos
        # and a mono that appeared at least 5 times, what percentage of the
        # top tetramers would we have found?
        d = sorted([(x, y) for x,y in counts.iteritems()], key=lambda z: z[1],
                   reverse=True)
        topGA_monos = set([x[0] for x in d[:20]])
        useful_monos = set([x for x, y in counts.iteritems() if y >= 4])
        assert topGA_monos.issubset(useful_monos)

        tot = 0
        for i, x in enumerate(toponepc):
            monos = x[1].split("_")[0].split("~")
            if monos[0] in useful_monos and monos[1] in useful_monos:
                if monos[0] in topGA_monos or monos[1] in topGA_monos:
                    tot += 1
        set_toponepc = set(x[1] for x in toponepc)
        ga_coverage = len(set_toponepc & res.pols)
        print "(%d, %d): %d found, %d out of %d would have been found" % (
            len(topGA_monos), len(useful_monos), ga_coverage, tot, len(toponepc))

        # Check that the best are in agreement
        
def setup_local_scan():
    length = 8
    mydata = {8:('octamer', 'octamers', 'afterga.txt', 70),
            6:('hexamer', 'hexamers', 'afterga_hexamers.txt', 50)}
    res = GA_Results(os.path.join(mydata[length][0], mydata[length][2]),
                     os.path.join(mydata[length][0], mydata[length][1]), length)
    pop = res.allpop

    data = {}

    for p in pop:
        data.update(dict(p))

    counts = collections.defaultdict(int)
    for x in data.keys():
        monos = x.split("_")[0].split("~")
        counts[monos[0]] += 1
        counts[monos[1]] += 1

    cutoff = mydata[length][3]
    best_monos = [x for x,y in counts.iteritems() if y >=cutoff]

    print "The total number of polymers generated was %d" % len(data)
    print "The number of monos that appear >= %d times is %d" % (cutoff, len(best_monos))

    output = list(counts.iteritems())
    output.sort(key=lambda x: x[1], reverse=True)
    for x, y in output:
        print x, y
    import pprint
    pprint.pprint([x[0] for x in output[:50]])

    pylab.hist(counts.values(), bins=range(0, 200, 2))
    pylab.show()

def plot_hex_and_oct_ga():
    """Create a Figure for the paper showing the Hexamer
    and Octamer GA and local search results"""

    folders = {6:'hexamer', 8:'octamer'}
    pcbm_lumo = -4.61
    admins = {6: os.path.join(folders[6], "hexamers"),
              8: os.path.join(folders[8], "octamers")}
    results = {6: os.path.join(folders[6], "afterga_hexamers.txt"),
               8: os.path.join(folders[8], "afterga_b.out")}

    zoomlevels = [0, 0, 2, 2]
    limits = {2: ((1.15, 2.4), (-7.1, -5.4)),
          1: ((0.5, 3.2), (-8.1, -4.8)),
          0: ((0.8, 6.5), (-10.5, -4.7)),
          }

    admin_dbs = dict([(x, Admin(y)) for x, y in admins.iteritems()])
    for x in admins.keys():
        print "The number of rows in the %d database is %d" % (x,
                                        len(admin_dbs[x].getalldata()))

    alldata = [None] * 4
    lengths = [6, 8, 6, 8]
    for i, length in enumerate([6, 8]):
        res = GA_Results(results[length], admins[length], length)
        data = set(res.pols)
        localpols = set([polname(x) for x in geneticAl.createAllCombinations(geneticAl.chosen_monomers[length],
                                                    length)])
        nlocalpols = localpols - data # Remove the existing ones
        print "%d: GA: %d, LS: %d, novel LS: %d" % (length, len(data), len(localpols),
                                                    len(nlocalpols))
        separate_localsearch_fromGA = False
        if separate_localsearch_fromGA:
            alldata[i] = data
            alldata[i + 2] = nlocalpols
        else:
            alldata[i] = data | nlocalpols
            alldata[i + 2] = alldata[i]

    print_info = True
    new3d = False # If you need to generate new 3D structures (slow)
    if print_info:
        # Print out the best hexamers + octamers
        mpl.rc('font', size=16)
        for N in [0, 1]:
            ax = pylab.subplot(121 + N)
            description = open(os.path.join(folders[lengths[N]], "top_description.txt"), "w")
            
            sc = ScoreCalculator()
            scores = []
            for name in alldata[N] | alldata[N + 2]:
                d = admin_dbs[lengths[N]].getdata(name)
                if d:
                    score = sc.getscore(d[1], cutoff=True)
                    scores.append( (name, score) )
            scores.sort(key=lambda x:x[1], reverse=True)

            print "Greater than 9%%", len([x for x in scores if x[1]>=9])
            print "Greater than 10%%", len([x for x in scores if x[1]>=10])
            print "Greater than 11%%", len([x for x in scores if x[1]>=11])
            
            a, b, c = pylab.hist([x[1] for x in scores], facecolor="gray", bins=np.arange(0.0, 12, 0.25))
            pylab.ylim(0, 250)
            pylab.title(folders[lengths[N]].capitalize())
            pylab.xlabel("Efficiency (%)")
            if N == 0:
                pylab.ylabel("Counts")
            ax.annotate("%d" % a[0], xy=(0.3, 235),  xycoords='data',
                xytext=(60, -40), textcoords='offset points',
                size=16,
                arrowprops=dict(arrowstyle="simple",
                                fc="0.6", ec="none",
                                patchB=c[0],
                                connectionstyle="arc3,rad=0.3"),
                )            
            
            if new3d:
                output = pybel.Outputfile("sdf", os.path.join(folders[lengths[N]], "top.sdf"), overwrite=True)
            for i in range(25):
                print >> description, "Top Number %d" % i
                print >> description, scores[i]
                mol = molname_to_mol(scores[i][0], lengths[N])
##                pdb.set_trace()
                print >> description, molname_to_repr(scores[i][0], lengths[N])
                json = admin_dbs[lengths[N]].getdata(scores[i][0])[1]
                homo, lumo, etens, etoscs = simplejson.loads(json)
                scale, trans, osc = besttrans(etens, etoscs, return_osc=True)                
                print >> description, "HOMO is %.3f" % homo
                print >> description, "LUMO is %.3f" % lumo
                print >> description, "Best trans is %.3f with osc strength %.3f" % (trans, osc)
                print >> description, "\n"
##                mol.draw(filename = os.path.join(folders[lengths[N]], "top%d.png" % i),
##                         show=False)
                if new3d:
                    globalopt(mol)
                    output.write(mol)
            if new3d:
                output.close()
        print pylab.gcf().get_size_inches()
        pylab.gcf().set_size_inches(9.5, 4)
        pylab.subplots_adjust(bottom=0.13, top=0.93, right=0.97, left=0.09)
        pylab.savefig("ForPaper_hex_oct_efficiency.png")

    plot = False
    if plot:
        figs = []
        
        CBs = []
        for N in range(4):
            figs.append(pylab.subplot(221 + N))

            x = []
            y = []
            fail = 0
            for name in alldata[N]:
                d = admin_dbs[lengths[N]].getdata(name)
                if d:
                    homo, trans = getHplusBG(d[1])
                    x.append(trans)
                    y.append(homo)
                else:
                    fail += 1
            pylab.plot(x, y, "g,")
            colorbar = False
            if N % 2 == 1:
                colorbar = True        
            CB = solar.plot(pcbm_lumo, colorbar=colorbar)
            if colorbar:
                CBs.append(CB)
                
            pylab.xlim(*limits[zoomlevels[N]][0]) 
            pylab.ylim(*limits[zoomlevels[N]][1])
            if N >= 2:
                pylab.xlabel("Lowest energy transition (eV)")
            if N % 2 == 0:
                pylab.ylabel("HOMO (eV)")
            if N < 2:
                pylab.title({6:"(a) Hexamers", 8:"(b) Octamers"}[lengths[N]])        

            print "The number of missing polymers out of %d is %d" % (len(alldata[N]), fail)

        pylab.subplots_adjust(right=0.95)
        lefts = [0.10, 0.52]
        bottom = [0.55, 0.07]
        bottom = [0.54, 0.08]
        width = [0.35, 0.35]
        height = [0.4, 0.4]
        figs[0].axes.set_position([lefts[0],bottom[0],width[0],height[0]])
        figs[1].axes.set_position([lefts[1],bottom[0],width[1],height[0]])
        figs[2].axes.set_position([lefts[0],bottom[1],width[0],height[1]])
        figs[3].axes.set_position([lefts[1],bottom[1],width[1],height[1]])
        CBs[0].ax.set_position([0.9,0.573,0.5,0.35])
        CBs[1].ax.set_position([0.9,0.086,0.5,0.35])
        for i in range(2):
            for t in CBs[i].ax.get_yticklabels():
                t.set_fontsize(9)
        pylab.savefig("Figure2.png")            

def plot_hex_and_oct_ga_heeger():
    """Create a Figure for the paper showing the Hexamer
    and Octamer GA and local search results"""

    supplementary = False
    
    folders = {6:'hexamer', 8:'octamer'}
    pcbm_lumo = -4.61
    admins = {6: os.path.join(folders[6], "hexamers"),
              8: os.path.join(folders[8], "octamers")}
    results = {6: os.path.join(folders[6], "afterga_hexamers.txt"),
               8: os.path.join(folders[8], "afterga_b.out")}

    if not supplementary:
        zoomlevels = [0, 0]
    else:
        zoomlevels = [1, 1]
    limits = {2: ((1.15, 2.4), (-7.1, -5.4)),
          1: ((6.5, 1.0), (-1.0, -6.2)),
          0: ((3.1, 1.0), (-3.3, -4.3)),
          }

    admin_dbs = dict([(x, Admin(y)) for x, y in admins.iteritems()])
    for x in admins.keys():
        print "The number of rows in the %d database is %d" % (x,
                                        len(admin_dbs[x].getalldata()))

    alldata = [None] * 4
    lengths = [6, 8, 6, 8]
    mpl.rc('font', size=16)
    for i, length in enumerate([6, 8]):
        res = GA_Results(results[length], admins[length], length)
        data = set(res.pols)
        localpols = set([polname(x) for x in geneticAl.createAllCombinations(geneticAl.chosen_monomers[length],
                                                    length)])
        nlocalpols = localpols - data # Remove the existing ones
        print "%d: GA: %d, LS: %d, novel LS: %d" % (length, len(data), len(localpols),
                                                    len(nlocalpols))
        separate_localsearch_fromGA = False
        if separate_localsearch_fromGA:
            alldata[i] = data
            alldata[i + 2] = nlocalpols
        else:
            alldata[i] = data | nlocalpols
            alldata[i + 2] = alldata[i]

    print_info = True
    new3d = False # If you need to generate new 3D structures (slow)
    if print_info:
        # Print out the best hexamers + octamers
        for N in [0, 1]:
            ax = pylab.subplot(121 + N)
            description = open(os.path.join(folders[lengths[N]], "top_description.txt"), "w")
            
            sc = ScoreCalculator()
            scores = []
            for name in alldata[N] | alldata[N + 2]:
                d = admin_dbs[lengths[N]].getdata(name)
                if d:
                    score = sc.getscore(d[1], cutoff=True)
                    scores.append( (name, score) )
            scores.sort(key=lambda x:x[1], reverse=True)

            print "Greater than 9%%", len([x for x in scores if x[1]>=9])
            print "Greater than 10%%", len([x for x in scores if x[1]>=10])
            print "Greater than 11%%", len([x for x in scores if x[1]>=11])
            
##            a, b, c = pylab.hist([x[1] for x in scores], facecolor="gray", bins=np.arange(0.0, 12, 0.25))
##            pylab.ylim(0, 250)
##            pylab.xlabel("%s Efficiency (%%)" % folders[lengths[N]].capitalize())
##            pylab.ylabel("Counts")
##            ax.annotate("Height is %d" % a[0], xy=(0.3, 235),  xycoords='data',
##                xytext=(60, -40), textcoords='offset points',
##                size=14,
##                arrowprops=dict(arrowstyle="simple",
##                                fc="0.6", ec="none",
##                                patchB=c[0],
##                                connectionstyle="arc3,rad=0.3"),
##                )            
##            
##            if new3d:
##                output = pybel.Outputfile("sdf", os.path.join(folders[lengths[N]], "top.sdf"), overwrite=True)
##            for i in range(25):
##                print >> description, "Top Number %d" % i
##                print >> description, scores[i]
##                mol = molname_to_mol(scores[i][0], lengths[N])
####                pdb.set_trace()
##                print >> description, molname_to_repr(scores[i][0], lengths[N])
##                json = admin_dbs[lengths[N]].getdata(scores[i][0])[1]
##                homo, lumo, etens, etoscs = simplejson.loads(json)
##                scale, trans, osc = besttrans(etens, etoscs, return_osc=True)                
##                print >> description, "HOMO is %.3f" % homo
##                print >> description, "LUMO is %.3f" % lumo
##                print >> description, "Best trans is %.3f with osc strength %.3f" % (trans, osc)
##                print >> description, "\n"
####                mol.draw(filename = os.path.join(folders[lengths[N]], "top%d.png" % i),
####                         show=False)
##                if new3d:
##                    globalopt(mol)
##                    output.write(mol)
##            if new3d:
##                output.close()
##        print pylab.gcf().get_size_inches()
##        pylab.gcf().set_size_inches(5, 8)
##        pylab.subplots_adjust(bottom=0.06, top=0.96, right=0.95)
##        pylab.savefig("ForPaper_hex_oct_efficiency.png")

    plot = True
    if plot:
        figs = []
        
        CBs = []
        for N in range(2):
            figs.append(pylab.subplot(121 + N))

            x = []
            y = []
            fail = 0
            for name in alldata[N]:
                d = admin_dbs[lengths[N]].getdata(name)
                if d:
                    homo, trans = getHplusBG(d[1])
                    x.append(trans)
                    y.append(homo+trans)
                else:
                    fail += 1
            pylab.plot(x, y, "g,")
            pylab.axhline(y = pcbm_lumo + 0.3, color="k", linewidth=1.0)            
            colorbar = False
            if N % 2 == 1:
                colorbar = True        
            CB = solar.plot_heeger(pcbm_lumo, colorbar=colorbar)
            if colorbar:
                CBs.append(CB)
                
            pylab.xlim(*limits[zoomlevels[N]][0]) 
            pylab.ylim(*limits[zoomlevels[N]][1])
            if N >= 0:
                pylab.xlabel("Transition energy (eV)")
            if N % 2 == 0:
                pylab.ylabel("Excited state (eV)")
            if N < 2:
                pylab.title({6:"Hexamers", 8:"Octamers"}[lengths[N]])        

            print "The number of missing polymers out of %d is %d" % (len(alldata[N]), fail)

        pylab.subplots_adjust(right=0.95)
        lefts = [0.10, 0.52]
        bottom = [0.11, 0.07]
        width = [0.35, 0.35]
        height = [0.8, 0.4]
        figs[0].axes.set_position([lefts[0],bottom[0],width[0],height[0]])
        figs[1].axes.set_position([lefts[1],bottom[0],width[1],height[0]])
##        figs[2].axes.set_position([lefts[0],bottom[1],width[0],height[1]])
##        figs[3].axes.set_position([lefts[1],bottom[1],width[1],height[1]])
        CBs[0].ax.set_position([0.91,0.15, 0.91, 0.75])
##        CBs[1].ax.set_position([0.9,0.086,0.5,0.35])
##        for i in range(1):
##            for t in CBs[i].ax.get_yticklabels():
##                t.set_fontsize(9)

        pylab.gcf().set_size_inches(9.5, 4)
        pylab.subplots_adjust(bottom=0.13, top=0.93, right=0.88, left=0.10)
        if not supplementary:
            pylab.savefig("Figure2.png")
        else:
            pylab.savefig("Supplem_Figure2.png")

def check_x_0_for_geoff():
    """Geoff: I was going through the analysis of the sequence
(i.e., I took the top 25 hexamers and octamers and set their
directions to 0,0).

It seems like the top results only come from the genetic algorithm,
and not the subsequent local search? I mention this, because in a
few cases, the 0,0 sequence had a higher efficiency than some of
the SMILES at the bottom of the list. I would assume the local search
would consider all directions of the top monomers, right?"""
    folders = {6:'hexamer', 8:'octamer'}
    admins = {6: os.path.join(folders[6], "hexamers"),
              8: os.path.join(folders[8], "octamers")}
    admin_dbs = dict([(x, Admin(y)) for x, y in admins.iteritems()])
    sc = ScoreCalculator()    
    for length in [6, 8]:
        print "*"*8, length, "*"*8
        for line in open(os.path.join(folders[length], "top_smiles.txt"), "r"):
            broken = line.split()
            eff = broken[-1]
            oldname = broken[0]
            polname = "_".join(oldname.split("_")[:2] + ["0"])

            print oldname, sc.getscore(admin_dbs[length].getdata(oldname)[1], cutoff=True)
            d = admin_dbs[length].getdata(polname)
            if d:
                score = sc.getscore(d[1], cutoff=True)
                print polname, score
            else:
                print polname, "Can't find it!"
            print



if __name__ == "__main__":
##    compare_objective_fns()
##    compare_chromos()
##    setup_local_scan()
    plot_hex_and_oct_ga_heeger()
##    plot_hex_and_oct_ga()
##    setup_conductivity()
##    check_x_0_for_geoff()
