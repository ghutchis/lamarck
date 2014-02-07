from cclib.parser import Gaussian

i = 1
f = open('../full-monomers.smi', 'r')
h = open('../homos-sort.txt', 'r')
e = open('../ex-sort.txt', 'r')
minhomo = 100.0
maxhomo = -100.0
minlumo = 100.0
maxlumo = -100.0
for line in f:
   homoLine = h.readline().split(' ')
   excitationLine = e.readline().split(' ')
   if len(homoLine) == 1:
      print i, line.rstrip()
      i += 1
      continue

   homo = float(homoLine[1])*27.21
   lumo = float(homo) + float(excitationLine[1])

   if homo < minhomo:
      minhomo = homo
   if homo > maxhomo:
      maxhomo = homo
   if lumo < minlumo:
      minlumo = lumo
   if lumo > maxlumo:
      maxlumo = lumo

   print i, line.rstrip(), homo, lumo, excitationLine[2].rstrip()
   i += 1

print "HOMO: ", minhomo, maxhomo
print "LUMO: ", minlumo, maxlumo
