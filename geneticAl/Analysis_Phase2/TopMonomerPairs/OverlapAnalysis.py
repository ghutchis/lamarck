#!/usr/bin/python
import sys
import os
import pandas as pd
import numpy as np
import pybel
import openbabel
import collections
from collections import OrderedDict

# Read in csv with top x% from different sized runs of data labeled SMILES_<number>
dataFile = sys.argv[1]
df = pd.read_csv(dataFile, sep=',', header=0)

x = (
    ('131', list(pd.Series(df['SMILES_131'].dropna()))),
    ('442', list(pd.Series(df['SMILES_442'].dropna()))),
    ('611', list(pd.Series(df['SMILES_611'].dropna()))),
    ('909', list(pd.Series(df['SMILES_909'].dropna()))),
    ('1235', list(pd.Series(df['SMILES_1235'].dropna()))),
    ('1759', list(pd.Series(df['SMILES_1759'].dropna())))
)

y = collections.OrderedDict(x)


intersection = {}
difference = {}
done = set()
created = []
for m in y:
    for n in y:
        k = "%s_%s" % (m, n)
        u = ''.join(sorted([m, n]))
        if (m is not n) and (u not in done):
            done.add(u)
            intersection[k] = map(str.strip, list(set(y[m]).intersection(y[n])))
            difference[k] = map(str.strip, list(set(y[m]).difference(y[n])))

            f_name_int = '%s_intersection' % k
            f_name_dif = '%s_difference' % k
            f_int = open('%s.smi' % f_name_int, 'w')
            f_dif = open('%s.smi' % f_name_dif, 'w')

            for i in intersection[k]:
                f_int.write("%s\n" % i)
            for i in difference[k]:
                f_dif.write("%s\n" % i)

            f_int.close()
            f_dif.close()

            created.append(f_name_int)
            created.append(f_name_dif)

# Makes a script for each of the combinations (intersections, differences, etc.) examined to quickly generate svg files
# by typing sh makesvg.sh
f_script = open('makesvg.sh', 'w')
f_script.write('#!/bin/sh\n')
for c in created:
    f_script.write('obabel %s.smi -O %s.svg --genalias -xA -xd -xb none -xC\n' % (c, c))
f_script.close()
