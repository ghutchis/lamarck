i = 0
f = open('octamer/top_smiles.txt', 'r')
for line in f:
   print line
   (data, efficiency) = line.rstrip().split(' ')
   (smiles, index1, index2) = data.split('_')

   nf = open('gaussian/' + str(i) + '.txt', 'w')
   nf.write(smiles + "_0_0") # change to 0,0 indexes -- most symmetric
   nf.close()
   i += 1
