Running the Genetic Algorithm
=============================

1. Create a similarity matrix from the monomers of interest which will be used in the geneticAl.py script to use
    the monomers of interest in the study. (These files are stored in the dims_and_tets folder):

    a. Create homopolymers of length 4 for each monomer of interest

       mkdir myhomopolymers
       python Homopolymer.py -f ../polysmiles.txt  -d myhomopolymers -l 4

    b. Run the PM6+ZINDO Gaussian calculations

       Go into myhomopolymers and run all of the .gjf Gaussian calculations

    c. Extract the results

        python ExtractData.py myhomopolymers

    d. Generate the similarity matrix

        python MakeSimMatrix.py myhomopolymers 4 mySimMatrix.json

        ** Make a note of the number of neighbours listed in the output **



2. Make a new folder for the output, say in geneticAl/output.
   
     python geneticAl.py --help

     python geneticAl.py -d output/myfirstrun --seed=0
                       --length=8   -N 64   
                       --matrix=mySimMatrix.json --nbrs=7
                       --objective=distance

The seed is the random seed, so that you can repeat the run exactly, or alternatively, try different runs.
The length is 8 (octamer) and there are 64 chromosomes.
The similiarty matrix used is stored in mySimMatrix.json (created as described in the README over in dims_and_tets),
    and the number of most similar neighbours to use is 7 (this value is given at the end of the creation of the similarity
    matrix, as described in Step 1.
The objective function is the distance.

Once this calculation has run, it can be repeated almost exactly by appending the --no-gaussian option,
which uses cached values instead of rerunning the Gaussian calculations. (It may not be exactly the same as some failed
calculations may have later succeeded.)
