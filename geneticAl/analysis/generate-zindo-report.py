#!/usr/bin/env python

import os
import sys
import json

def exitwitherror():
    sys.exit(__doc__)

minhomo = -10.72
midhomo = -8.537
maxhomo = -6.94

def homocolor(energy):
   # Retrieved from monomer DFT calculations
   # "mid" is reference -- polythiophene
   if energy < midhomo:
      # red to white
      range = midhomo - minhomo
      distance = energy - minhomo
      index = int(distance/range*255)
      color = '%x' % index
      if len(color) == 1:
          color = '0' + color
      return 'ff%s%s' % (color, color)
   else:
      # white to blue
      range = maxhomo - midhomo
      distance = energy - midhomo
      index = 255 - int(distance / range * 255)
      color = '%x' % index
      if len(color) == 1:
          color = '0' + color
      return '%s%sff' % (color, color)

minlumo =  -6.60
midlumo =  -2.03
maxlumo =  -0.393
def lumocolor(energy):
   if energy < midlumo:
      # yellow to white
      range = midlumo - minlumo
      distance = energy - minlumo
      index = int(distance/range*255)
      color = '%x' % index
      if len(color) == 1:
          color = '0' + color
      return 'ffff%s' % color
   else:
      # white to green
      range = maxlumo - midlumo
      distance = energy - midlumo
      index = 255 - int(distance / range * 255)
      color = '%x' % index
      if len(color) == 1:
          color = '0' + color
      return '%sff%s' % (color, color)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        exitwitherror()

    print "<html><head><title>Analysis Results</title></head><body>"

    print "<h2>Table S3:</h2>"
    print "<p>Top octamer results, including predicted efficiencies, energetics of component momomers, and sequence effects.</p>"
    print "<table border=0>"
    print "<tr><th>Global Color Scales</th><th>Minimum</th><th>Mid-Point</th><th>Maximum</th></tr>"
    print "<tr><th>HOMO Scale (in eV)</th>"
    print "<td bgcolor=\"" + homocolor(minhomo) + "\">" + "%6.2f" % minhomo + "</td>"
    print "<td bgcolor=\"" + homocolor(midhomo) + "\">" + "%6.2f" % midhomo + "</td>"
    print "<td bgcolor=\"" + homocolor(maxhomo) + "\">" + "%6.2f" % maxhomo + "</td>"
    print "</tr>"
    print "<tr><th>LUMO Scale (in eV)</th>"
    print "<td bgcolor=\"" + lumocolor(minlumo) + "\">" + "%6.2f" % minlumo + "</td>"
    print "<td bgcolor=\"" + lumocolor(midlumo) + "\">" + "%6.2f" % midlumo + "</td>"
    print "<td bgcolor=\"" + lumocolor(maxlumo) + "\">" + "%6.2f" % maxlumo + "</td>"
    print "</tr>"
    print "</table><p>&nbsp;</p>"

    print "<table border=0>"

    filename = sys.argv[1]
    oligolen = sys.argv[2]

    combs = json.load(open("directions.json", 'r'))
    # combs['8'][#][#]
    # (di)(qu)(di)(qu)(di)(qu)(di)(qu)

    monodatafile = open("monomer-zindo.txt", 'r')
    homos = {}
    lumos = {}
    for line in monodatafile:
        dataList = line.split(' ')
        if len(dataList) < 5:
            continue
        # (number, smiles, homo, lumo, oscstrs)
        smiles = dataList[1]
        homos[smiles] = float(dataList[2])
        lumos[smiles] = float(dataList[3])

    bestdimers = open("best-dimers.smi", 'r')
    dimerindex = {}
    index = 0
    for line in bestdimers:
        smiles = line.rstrip()
        dimerindex[smiles] = index
        index += 1

    f = open(filename, 'r')
    print "<tr><th>Index</th><th>Efficiency</th><th>HOMO (eV)</th><th>LUMO (eV)</th><th colspan=\"" + str(oligolen) + "\">HOMO Barcode</th><th colspan=\"" + str(oligolen) + "\">LUMO Barcode</th></tr>"
    entry = 1
    for line in f:
       print "<tr>"
       (data, efficiency, homo, excit) = line.rstrip().split(' ')
       lumo = float(homo) + float(excit)
       (smiles, index1, index2) = data.split('_')
       (mono1, mono2) = smiles.split('~')
       pattern = combs[oligolen][int(index1)][int(index2)]
       segments = pattern[1:-1].split(')(')

       try:
           index = dimerindex[smiles]
       except:
           index = 0
#       print "<td>" + '%5.2f' % float(efficiency) + "</td><td>" + "<img src='../../donoracceptor/dimerdepictions/%d.png' height=128>" % index + "</td>"
       print "<td>%d</td>" % entry + '<td>%5.2f</td><td>%5.2f</td><td>%5.2f</td>' % (float(efficiency), float(homo), float(lumo))
       # code up homo colors
       for segment in segments:
           if segment == 'uq' or segment == 'qu':
               print "<td bgcolor=\"" + homocolor(homos[mono1]) + "\">" + "%6.2f" % homos[mono1] + "</td>"
           else:
               print "<td bgcolor=\"" + homocolor(homos[mono2]) + "\">" + "%6.2f" % homos[mono2] + "</td>"
       for segment in segments:
           if segment == 'uq' or segment == 'qu':
               print "<td bgcolor=\"" + lumocolor(lumos[mono1]) + "\">" + "%6.2f" % lumos[mono1] + "</td>"
           else:
               print "<td bgcolor=\"" + lumocolor(lumos[mono2]) + "\">" + "%6.2f" % lumos[mono2] + "</td>"
       print "</tr>"
       entry += 1
    f.close()

    print "</table>"
    print "</body></html>"
