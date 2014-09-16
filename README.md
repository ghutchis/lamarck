Running the Genetic Algorithm
=============================
1. Make a library of monomers to include in the screening:
    
      a. The library should be a txt file with one monomer per line with the beginning and ending letters of the
        code being the polymerization sites so that when the string is placed beside itself 4-8 times (tetramers- octomers) the SMILES do not produce errors. An example list is shown here:
        c(s1)c(S(=O)(=O)C=C2)c2c1
        c(s1)c(OCCCO2)c2c1
        c(o1)nnc1
        c(c(nsn1)c12)ccc2
        c(s1)c(SC=CS2)c2c1
            etc.....
    
    b. Generate gaussian input files for homotetramers (using input.py file and specify calculation type in the       header) and run gaussian calcuations for each monomer.

    c. Find all cations to determine the probable polymerization sites.
        1) TODO: Is this done with the g09catreorg.sh script on the server?

    d. Generate complete list on monomers with all possible polymerization sites included.


2. Create a similarity matrix from the monomers of interest which will be used in the geneticAl.py script to use
    the monomers of interest in the study. This should be completed for each length of interest. The directions
    describe the process for tetramers, but can be changed for dimers, hexamers, octomers, etc. (These files are stored in the dims_and_tets folder):

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



3. Run genetic algorithm:

    a. Make a new folder for the output, say in geneticAl/output.
   
    b. The syntax to run the genetic algorithm is:
        python geneticAl.py -d output/<run name> --seed=0 --length=4   -N 64 --matrix=mySimMatrix.json --nbrs=7 --objective=distance

        1) The seed is the random seed, so that you can repeat the run exactly, or alternatively, try different runs.
        2) The length is 4 (tetramer) and there are 64 chromosomes.
        3) The similiarty matrix used is stored in mySimMatrix.json (created as described in previous step and
            corresponding to the length chosen.
        4) Number of most similar neighbors to use is 7 (this value is given at the end of the creation of the
            similarity matrix, as described in previous step.
        5) The objective function is the distance.

    c. If help is needed, use: python geneticAl.py --help

    d. Once this calculation has run, it can be repeated almost exactly by appending the --no-gaussian option,
        which uses cached values instead of rerunning the Gaussian calculations. (It may not be exactly the same as
        some failed calculations may have later succeeded.)

    e. To submit a job to the Hutchison clusher:
        1. Adjust the keywords in the runGeneticAl.sh file
        2. Type: qsub -N <job name>  /Volumes/Users/ilanakanal/screeningproject/geneticAl/RunGeneticAl.sh
    f. To preselect the initial population, add a set of monomers from which the initial set should be chosen. In addition,
        CHANGE:  monos = [random.choice(self.monomers) for j in range(2)]
        TO:      monos = [random.choice(selected_initial_population) for j in range(2)]
        and CHANGE:  def initallpop(self, chosen_monos=None):
            TO:      def initallpop(self, chosen_monos=selected_initial_population):

        When generating the list of monomers to be included in the selected_initial_population list, care must be taken
         to make sure each monomer is in the similarity matrix which is used in the run (otherwise the next generation
          will not be able to be generated and the program will fail.

            Example: selected_initial_population = ['c(s1)c(cc(C#N)c(C#N)c2)c2c1','c(s1)c(C(=O)OC2(=O))c2c1']

4. Processing data:
    a. Open the desired .db file which can be found in the geneticAl/output folder using a SQL database browser.
        Save file as a (csv) text file. Delete the top row (header row). SAVE this text file in:
         geneticAl/Analysis_Phase2