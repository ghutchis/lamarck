from cclib.parser import Gaussian
from Efficiency import *
import json
import logging

efficient = Efficiency()

h = open('hexamer/hexamersDB.txt', 'r')
homos = {}
exState = {}
exStr = {}
for line in h:
   data = line.split('"')
   smiles = data[1].split('_')[0]
   eData = json.loads(data[5])
   homos[smiles] = float(eData[0])

   eTrans = eData[2]
   eStr = eData[3]

   bestExcitationIdx = 0
   bestExcitationStr = 0.0
   j = 0
   for excitation in eStr:
      if excitation > bestExcitationStr:
         bestExcitationStr = excitation
         bestExcitationIdx = j
      j += 1

   # OK, now bestExcitationIdx has the best index
   exState[smiles] = homos[smiles] + float(eTrans[bestExcitationIdx])
   exStr[smiles] = bestExcitationStr

i = 0
t = open('hexamer/top_smiles.txt', 'r')

for line in t:
   smiles = line.split(' ')[0].split('_')[0]
   eff = line.split(' ')[1].rstrip()

   myfile = Gaussian("gaussian/hex-00/%d.g09" % i)
   myfile.logger.setLevel(logging.ERROR)
   data = myfile.parse()
   homoE = data.moenergies[0][data.homos[0]]
   excitation = data.etenergies[0] * 0.00012398
   lumoE = homoE + excitation
   
   eff2 = efficient.efficiency(homoE, excitation, -4.61)

   print smiles, eff, homos[smiles], exState[smiles], exStr[smiles], eff2, homoE, lumoE, data.etoscs[0]
   i += 1
