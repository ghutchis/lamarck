#!/usr/bin/env python

import sys
from Efficiency import *

efficient = Efficiency()

if (len(sys.argv) == 2):
    for line in open(sys.argv[1]):
        (homo, bandgap) = line.split()
        print efficient.zindoEff(float(homo), float(bandgap))
elif (len(sys.argv) > 3):
    print "B3LYP Efficiency: ", efficient.b3lypEff(float(sys.argv[1]), float(sys.argv[2]))
else:
    print "ZINDO Efficiency: ", efficient.zindoEff(float(sys.argv[1]), float(sys.argv[2]))
