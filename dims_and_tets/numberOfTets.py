import sys

N = 66
sim = [[x, x] for x in range(N)]
unsim = [[x, x+1] for x in range(N)]
allmonos = sim + unsim + [[x[1], x[0]] for x in unsim]

def polname(*vals):
    pol = list(vals)
    rev = [[x[1], x[0]] for x in pol]
    rev.reverse()
    if pol < rev:
##        return "-".join(str(x) for x in pol)
        return str(pol)
    else:
        return str(rev)
##        return "-".join(str(x) for x in rev)

def dimers():
    # ===== Dimers =====
    out = open("tmp.txt", "w")
    for x in allmonos:
        for y in allmonos:
            print >> out, polname(x, y)
    out.close()
    # Run out through sort and uniq to find the number of dimers
    # Verify equal to 19701
    # /usr/bin/sort -i tmp.txt | wc -l

def tetramers():
    # ===== Tetramers =====
    out = open("tmp_b.txt", "w")
    for i, x in enumerate(allmonos):
        print i, 
        for y in allmonos:
            for z in allmonos:
                for za in allmonos:
                    print >> out, polname(x, y, z, za)
    out.close()

if __name__ == "__main__":
    dimers()
    tetramers()
