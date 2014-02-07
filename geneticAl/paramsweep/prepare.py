template=  "python ../geneticAl.py --database=/ichec/home/users/noboyle/Work/solar/dims_and_tets/%s --no-gaussian --length=%d -N %d --nbrs=%d --matrix=%d --objective=%s --seed=%d > output/%s_CHR%d_NBR%d_MAT%d_%s_SEED%d.txt"

dbs = {2: "alldimers", 4: "alltetramers"}

lines = set()
for length in [2, 4]:
    nchromo = 64
    nbrs = 7
    matrix = 4
    seed = 0
    for objfn in ['eff', 'distance']:
        for matrix in [2, 4]:
            for seed in range(10):
                lines.add(template % (dbs[length], length, nchromo, nbrs, matrix,
                      objfn, seed, dbs[length], nchromo, nbrs, matrix, objfn, seed))
        seed = 0
        matrix = 4
        for nchromo in [8, 16, 32, 64, 128, 256]:
            lines.add(template % (dbs[length], length, nchromo, nbrs, matrix,
                  objfn, seed, dbs[length], nchromo, nbrs, matrix, objfn, seed))
        nchromo = 64
        for nbrs in [3, 7, 11, 15, 31, 66]:
            lines.add(template % (dbs[length], length, nchromo, nbrs, matrix,
                  objfn, seed, dbs[length], nchromo, nbrs, matrix, objfn, seed))
        nbrs = 7


output = open("tasks", "w")
for line in lines:
    print >> output, line
output.close()
