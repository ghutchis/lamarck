__author__ = 'ilanakanal'
#!/usr/bin/python
import sys
from collections import deque
import json
import json as simplejson

# Argument is the text file with the monomers which should be combined to make all possible symmetry combinations for each sequence.
# Files which are to be used are SequenceMonomers.txt, SequenceMonomers2.txt, SequenceMonomers3.txt and SequenceMonomers4.txt for
# both assymmetric, first monomer is symmetric, second monomer is symmetric and both are symmetric. To use the file, change the monomer
# length to the desired length and then copy the output to the directions.json file

#TODO In order to get the pairwise combinations to form given combinations, there are repeats (ie: AAAAAA, etc). Does this need to be changed so that the likelihood of picking a given combination is the same?


datafile = open(sys.argv[1])

for line in datafile.readlines():
    line = line.rstrip()
    #name = "line" + line.split()[1] + line.split()[3]
    donor = line.split()[0]
    acceptor = line.split()[1]
    length = 6
    mon1 = 'D'
    mon2 = 'A'

    Sequences = []
    data = {}
    finalList = ''
    list = []
    final = ''
    SeqList = []

    def combinations(iterable, r):
        pool = tuple(iterable)
        n = len(pool)
        if r > n:
            return
        indices = range(r)
        yield tuple(pool[i] for i in indices)
        while True:
            for i in reversed(range(r)):
                if indices[i] != i + n - r:
                    break
            else:
                return
            indices[i] += 1
            for j in range(i + 1, r):
                indices[j] = indices[j - 1] + 1
            yield tuple(pool[i] for i in indices)

    def tf_gen(N):
        stack = deque()
        stack.append([])
        while len(stack):
            elem = stack.pop()
            if len(elem) == N - 1:
                yield elem + [False]
                yield elem + [True]
            else:
                stack.append(elem + [True])
                stack.append(elem + [False])

    SeqList = tf_gen(length)

    for length in [2, 4, 6, 8]:
       #print "\n****", length, "****"

        # Does not filter out reverse or 50:50 mixes
        for Seq in SeqList:
            SeqString0 = ''
            for item in Seq:
                if (item):
                    SeqString0 += mon1
                else:
                    SeqString0 += mon2
            data[SeqString0] = 1

            SeqString = SeqString0.replace(mon1, donor).replace(mon2, acceptor)
            print ""+'"%s"' % SeqString+ ',' # Displays possible donor/ acceptor combinations

            #with open("directions_sequences2.json", "w") as f:
#                simplejson.dump(""+SeqString+  ',', f)