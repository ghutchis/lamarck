"""Create directions.pickle

This file contains all of the possible ways to combine a particular
monomer with another to create polymers of different length,
taking into account their symmetries and the
requirement that each polymer is built up from a base dimer.
"""

import pdb
import pybel
from collections import deque
import json as simplejson

reverse = {"C=N":"N=C", "SO":"OS", 'C=C':'C=C', 'SS':'SS'}

def tf_gen(N):
    """Generator of all possible combinations of True/False"""
    stack = deque()
    stack.append([])
    while len(stack):
        elem = stack.pop()
        if len(elem) == N-1:
            yield elem + [False]            
            yield elem + [True]
        else:
            stack.append( elem+[True] )
            stack.append( elem+[False] )

def join(smiles, directions):
    ans = []
    for x,y in zip(smiles, directions):
        if y:
            ans.append(x)
        else:
            ans.append(reverse[x])
    smile = "".join(ans)
    return smile

def main():
    differentmonos = [["C=N", "SO"],["C=C", "SO"],["C=N", "SS"],["C=C", "SS"],
                      ["C=N", "C=N"], ["C=C", "C=C"]]

    data = {}
    for length in [2, 4, 6, 8]:
        print "\n****", length, "****"
    
        allresults = []
        for a in differentmonos:
            results = []
            already = set()
            for d in tf_gen(length/2 + 3): # 4 bits for 2mers, 5 for 4mers
                # Set the basic unit
                b = [a[0], a[1]]
                if d[0]:
                    b = [a[1], a[0]]

                smiles = []
                directions = []
                for x in d[3:]:
                    if x:
                        directions.extend([d[1], d[2]])
                        smiles.extend([b[0], b[1]])
                    else:
                        directions.extend([not d[2], not d[1]])
                        smiles.extend([b[1], b[0]])
                newsmiles = join(smiles, directions)
                can = pybel.readstring("smi", newsmiles).write("can").split()[0]
                if can not in already:
                    result = "".join(newsmiles).replace(
                        a[0], "qu").replace(reverse[a[0]], "uq").replace(a[1], "di").replace(
                        reverse[a[1]], "id")
                    print smiles, directions, newsmiles, result
                    template = "(%s)" % ")(".join(result[i:i+2] for i in range(0, len(result), 2))
                    results.append(template)
                    already.add(can)
            allresults.append(results)
            print

        data[length] = allresults
        print
    for length in [2, 4, 6, 8]:
        assert len(data[length][-1][0]) / 4 == length
        lengths = [len(x) for x in data[length]]
        print lengths
        print (lengths[0]*(66*65 / 2) + lengths[1]*66*66 +
               lengths[3]*(66*65/2) + lengths[4]*66 + lengths[5]*66)

    with open("directions.json", "w") as f:
        simplejson.dump(data, f)

if __name__ == "__main__":
    main()
