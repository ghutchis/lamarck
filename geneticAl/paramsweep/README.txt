What parameters are best for the genetic algorithm? Does it make any difference whether we only mutate to the 7 nearest neighbors, or whether we randomly choose where to mutate to.

Let's find out.

Since we know the complete space of dimers and tetramers, we can investigate different parameter values and see how they compare. Later, for hexamers and octamers, we will simply use our estimate of the best parameters.

python geneticAl.py --help gives:
  -s SEED, --seed=SEED  set the random seed
  -d DATABASE, --database=DATABASE
                        set the database name
  -t, --test            run the test
  --no-gaussian         don't run any Gaussian jobs
  -r, --random          chose the next generation randomly
  -l LENGTH, --length=LENGTH
                        set the polymer length
  -N SIZE               set the number of chromosomes
  -R NBRS, --nbrs=NBRS  set the number of nearest neighbours
  --matrix=MATRIX_ID    which similarity matrix to use (2 or 4)
  --objective=OBJ_FN    which objective function (eff or distance)

(1) By altering the seed from 1 to 10 we can repeat the same run several times with different initial populations.
(2) By altering the database, we can focus on dimers or tetramers. Remember to set the appropriate length value (i.e. 2 for dimers) though or the results are meaningless.
(3) We will alter the number of chromosomes: 8, 16, 32, 64, 128 or 256.
(4) We will look at using 4, 7, 11, 15, 19 or 70 neighbors.
(5) Should we use the similarity matrix based on the homodimers or homotetramers (SIM2 or SIM4)?
(6) Which objective function: efficiency above all else, or distance to the known sweet spot?
