Running the Genetic Algorithm
=============================

1. Make a new folder for the output, say in geneticAl/output.
   
     python geneticAl.py --help

     python geneticAl.py -d output/myfirstrun --seed=0
                       --length=8   -N 64   
                       --matrix=mostsim4.json --nbrs=7
                       --objective=distance

The seed is the random seed, so that you can repeat the run exactly, or alternatively, try different runs.
The length is 8 (octamer) and there are 64 chromosomes.
The similiarty matrix used is stored in mostsim4.json (created as described in the README over in dims_and_tets), and the number of most similar neighbours to use is 7.
The objective function is the distance.

Once this calculation has run, it can be repeated almost exactly by appending the --no-gaussian option, which uses cached values instead of rerunning the Gaussian calculations. (It may not be exactly the same as some failed calculations may have later succeeded.)
