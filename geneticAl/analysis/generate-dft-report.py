#!/usr/bin/env python

import os
import sys
import json

def exitwitherror():
    sys.exit(__doc__)


def homocolor(energy):
   # Retrieved from monomer DFT calculations
   # "mid" is reference -- polythiophene
   minhomo =  -8.58628
   midhomo =  -5.66213
   maxhomo =  -3.84796
   if energy < midhomo:
      # red to white
      range = midhomo - minhomo
      distance = energy - minhomo
      index = int(distance/range*255)
      return 'ff%x%x' % (index, index)
   else:
      # white to blue
      range = maxhomo - midhomo
      distance = energy - midhomo
      index = 255 - int(distance / range * 255)
      return '%x%xff' % (index, index)

def lumocolor(energy):
   minlumo =  -3.93449
   midlumo =  -0.24380
   maxlumo =   1.38642
   if energy < midlumo:
      # yellow to white
      range = midlumo - minlumo
      distance = energy - minlumo
      index = int(distance/range*255)
      return 'ffff%x' % index
   else:
      # white to green
      range = maxlumo - midlumo
      distance = energy - midlumo
      index = 255 - int(distance / range * 255)
      return '%xff%x' % (index, index)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        exitwitherror()

    print "<html><head><title>Analysis Results</title></head><body>"
    print "<table border=0>"

    filename = sys.argv[1]
    oligolen = sys.argv[2]

    combs = json.load(open("directions.json", 'r'))
    # combs['8'][#][#]
    # (di)(qu)(di)(qu)(di)(qu)(di)(qu)

    monodatafile = open("monomer-dft.txt", 'r')
    homos = {}
    lumos = {}
    monomers = {}
    for line in monodatafile:
        (number, smiles, homo, lumo) = line.split(' ')
        homos[smiles] = float(homo)
        lumos[smiles] = float(lumo)
        monomers[smiles] = "l" + number + ".g03.gz"

    f = open(filename, 'r')
    print "<tr><th>Efficiency</th><th colspan=\"" + str(oligolen) + "\">HOMO Barcode</th><th colspan=\"" + str(oligolen) + "\">LUMO Barcode</th></tr>"
    for line in f:
       print "<tr>"
       (data, efficiency) = line.rstrip().split(' ')
       (smiles, index1, index2) = data.split('_')
       (mono1, mono2) = smiles.split('~')
       pattern = combs[oligolen][int(index1)][int(index2)]
       segments = pattern[1:-1].split(')(')

       print "<td>" + '%5.2f' % float(efficiency) + "</td>"
       # code up homo colors
       for segment in segments:
           if segment == 'uq' or segment == 'qu':
               print "<td bgcolor=\"" + homocolor(homos[mono1]) + "\"><a href=\"" + monomers[mono1] + "\">&nbsp;</a></td>"
           else:
               print "<td bgcolor=\"" + homocolor(homos[mono2]) + "\">&nbsp;</td>"
       for segment in segments:
           if segment == 'uq' or segment == 'qu':
               print "<td bgcolor=\"" + lumocolor(lumos[mono1]) + "\">&nbsp;</td>"
           else:
               print "<td bgcolor=\"" + lumocolor(lumos[mono2]) + "\">&nbsp;</td>"
       print "</tr>"
    f.close()

    print "</table>"
    print "</body></html>"
