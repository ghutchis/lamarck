import os
import pdb
import math
import pylab
import numpy as np
try:
    import json
except ImportError:
    import simplejson as json
    
import collections

import solar
    
from Efficiency import *
from Utils import *
from Admin import *

relpath = os.sep.join(__file__.split(os.sep)[:-1])

def within(distance, (x1, y1), (x2, y2)):
    d = math.sqrt((x1-x2)**2 +(y1-y2)**2)
    return d < distance

def closetoline(distance, (x0, y0), (x1, y1), (x2, y2)):
    d = abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)) / math.sqrt(
                (x2-x1)**2 + (y2-y1)**2)
    return d < distance

def zoomin():
    """Plot the full set of dimers and tetramers along with a
    zoomed-in version"""
    
    sc = ScoreCalculator()
    pcbm_lumo = -4.61

    dbname = {2: 'alldimers', 4: 'alltetramers'}
    lengths = [2, 4, 2, 4]
    zoomlevels = [0, 0, 2, 2]
    limits = {2: ((1.15, 2.4), (-7.1, -5.4)),
          1: ((0.5, 3.2), (-8.1, -4.8)),
          0: ((0.8, 6.5), (-10.5, -4.7)),
          }
    alldata = {2: Admin("alldimers").getalldata(),
               4: Admin("alltetramers").getalldata()}

    recalc = True # The first time you run this you need to use recalc = True
    alleffs = {2:[], 4:[]}
    if not recalc:
        alleffs = json.load(open("tmp.json", "r"))
        alleffs[2] = alleffs['2']
        alleffs[4] = alleffs['4']
        del alleffs['2']
        del alleffs['4']

    for length in [2, 4]:
        ax = pylab.subplot(210 + length/2)
        if recalc:
            for d in alldata[length]:
                homo, trans = getHplusBG(d[2])
                e = sc.getscore(d[2])
                alleffs[length].append(e)

        a, b, c = pylab.hist(alleffs[length], facecolor="gray", bins=np.arange(0.0, 12, 0.25))
        pylab.ylim(0, 2750)        
        pylab.xlabel("%s Efficiency (%%)" % dbname[length][3:-1].capitalize())
        pylab.ylabel("Counts")

        ax.annotate("Height is %d" % a[0], xy=(0.3, 2580),  xycoords='data',
                xytext=(60, -40), textcoords='offset points',
                size=14,
                arrowprops=dict(arrowstyle="simple",
                                fc="0.6", ec="none",
                                patchB=c[0],
                                connectionstyle="arc3,rad=0.3"),
                )
        # Move tick marks outside y-axis
        lines = ax.get_yticklines()
        labels = ax.get_yticklabels()
        for line in lines:
            line.set_marker(pylab.matplotlib.lines.TICKLEFT)
        for label in labels:
            label.set_x(-0.02)        

    if recalc:
        json.dump(alleffs, open("tmp.json", "w"))
    
    print pylab.gcf().get_size_inches()
    pylab.gcf().set_size_inches(5, 8)
    pylab.subplots_adjust(bottom=0.06, top=0.96, right=0.95, left=0.16)
##    pylab.show()
    pylab.savefig(os.path.join("pictures", "ForPaper_hist_efficiency.png"))

    figs = []
    CBs = []
    for N in range(4):
        figs.append(pylab.subplot(221 + N))
        length = [2, 4, 2, 4][N]
        db = dbname[length]
        x = []
        y = []
        selected = []
        for d, e in zip(alldata[length], alleffs[length]):
            homo, trans = getHplusBG(d[2])
            highlight = False
            if within(2, (homo, trans), (-5.75, 1.45)):
##                if e > 6.0:
                if e > 8.0:
                    highlight = True
                elif (e == 0 and closetoline(
                        0.1, (trans, homo), (2.0, -6.2), (1.15, -5.7))
                      and sc.getscore(d[2], cutoff=False) > 8.0
                    ):
                    highlight = True
            if highlight:
                selected.append( (d[0], molname_to_smile(d[0], length), trans, homo, e))
            else:               
                y.append(homo)
                x.append(trans)
        output = open("best_%s.txt" % db, "w")
        for s in selected:
            print >> output, "\t".join([str(noel) for noel in s])
        output.close()

        zoomlevel = zoomlevels[N]
        dot = ","
        if zoomlevel == 2:
            dot = "."
        pylab.plot(x, y, "g" + dot)
        pylab.plot([q[2] for q in selected],
                   [q[3] for q in selected], "r" + dot)
        colorbar = False
        if N % 2 == 1:
            colorbar = True
        CB = solar.plot(pcbm_lumo, colorbar=colorbar)
        if colorbar:
            CBs.append(CB)

        pylab.xlim(*limits[zoomlevel][0])
        pylab.ylim(*limits[zoomlevel][1])
        if N >= 2:
            pylab.xlabel("Lowest energy transition (eV)")
        if N % 2 == 0:
            pylab.ylabel("HOMO (eV)")
        if N < 2:
            pylab.title({2:"(a) Dimers", 4:"(b) Tetramers"}[length])
    pylab.subplots_adjust(right=0.95)
    lefts = [0.10, 0.52]
    bottom = [0.55, 0.07]
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
            
    pylab.savefig(os.path.join("pictures", "Figure1_dims_and_tets.png"))
##    pylab.show()

def hist():
    """Draw histograms of the distribution of HOMOs and transition Es
    for the dimers vs the tetramers"""
    
    eff = Efficiency()
    pcbm_lumo = -4.61
    homos = [[], []]
    trans = [[], []]
    for i in range(2):
        db = ['alldimers', 'alltetramers'][i]
        admin = Admin(db)        
        x = []
        y = []
        selected = []
        for d in admin.getalldata():
            homo, tran = getHplusBG(d[2])
            homos[i].append(homo)
            trans[i].append(tran)
        print np.mean(homos[i]), np.mean(trans[i])
    delta = 0.2
    bins = [np.arange(-10, -5.5, delta), np.arange(1, 7, delta)]

    titles = ['HOMO (eV)', 'Lowest energy significant transition (eV)']
    for i in range(2):
        data = [homos, trans][i]
        n, bns = np.histogram(data[1], bins[i])
        c = pylab.bar(bns[:-1], n/float(n.sum()), delta * .4, color="gray")

        for rect in c:
            rect.set_x(rect.get_x() + delta/2. * 0.8)
        n, bns = np.histogram(data[0], bins[i])
        f = pylab.bar(bns[:-1], n/float(n.sum()), delta * .4, color="k")

        pylab.legend([f[0], c[0]], ["Dimers", "Tetramers"])
        pylab.ylabel("Fraction")
        pylab.xlabel(titles[i])
        pylab.savefig(os.path.join("pictures", "Figure1_hist_%d.png" % i))
        pylab.clf()

def most_eff_tetramers():
    data = Admin("alltetramers").getalldata()
    sc = ScoreCalculator()
    res = []
    for d in data:
        smiles = molname_to_smile(d[0], 4)
        score = sc.getscore(d[2])
        res.append((score, d[0], smiles))
    res.sort(reverse=True)
    print "\n".join(str(x) for x in res[0:20])

def best_monomers():
    """How many monomers were involved in the top N solutions?"""
    db = {2: "alldimers", 4: os.path.join(relpath, "alltetramers")}
    length = 4
    data = Admin(db[length]).getalldata()
    sc = ScoreCalculator()
    res = []
    for d in data:
        smiles = molname_to_smile(d[0], length)
        dist = sc.getdistance(d[2])
        res.append((dist, d[0], smiles))
    res.sort()

    ans = []
    total = set()
    ic = 0
    cutoffs = [16, 33, 66, 10000]
    counts = collections.defaultdict(int)
    onepc = len(data) / 100
    for i, x in enumerate(res):
        monos = x[1].split("_")[0].split("~")
        counts[monos[0]] += 1
        counts[monos[1]] += 1        
        total.update(monos)
        ans.append(len(total))
        if ans[-1] > cutoffs[ic]:
            print "%d monos cover the top %d" %  (cutoffs[ic], i)
            ic += 1

        if i == onepc:
            print "Top 1%% is %d and involves %d monos" % (len(data)/100, ans[-1])
            counts_onepc = counts.copy()

    bestmonos = set([x for x, y in counts_onepc.iteritems() if y > 40])
    tot = 0
    for i, x in enumerate(res[:onepc]):
        monos = x[1].split("_")[0].split("~")
        if monos[0] in bestmonos or monos[1] in bestmonos:
            tot += 1
    print "%d out of %d have one of the best %d monos" % (tot, onepc, len(bestmonos))
    
##    a = pylab.hist(counts_onepc.values(), bins=range(0, 100, 2))
##    pylab.show()

##     bestmonos = set([u'c(s1)c(NC=N2)c2c1', u'c(s1)cc(c12)Cc(c2s3)cc3', u'c(s1)c(nccn2)c2c1',
##                      u'c(s1)c(ncnc2)c2c1', u'c1oc(c2)c(c1)oc2',
##                      u'c(s1)cc(c12)[nH]c(c2s3)cc3', u'c(s1)c(NC=C2)c2c1'])
    return bestmonos, res[:onepc]
  

if __name__ == "__main__":
    zoomin()
##    best_monomers()
##    hist()

