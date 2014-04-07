Creating the Similarity Matrix
==============================

Summary:

1. Create homopolymers of length 4 for each monomer of interest
2. Run the PM6+ZINDO Gaussian calculations
3. Extract the results
4. Generate the similarity matrix

In full:

1. Create homopolymers of length 4 for each monomer of interest

   mkdir myhomopolymers
   python Homopolymer.py -f ../polysmiles.txt  -d myhomopolymers -l 4

2. Run the PM6+ZINDO Gaussian calculations

   Go into myhomopolymers and run all of the .gjf Gaussian calculations

3. Extract the results

   python ExtractData.py myhomopolymers

4. Generate the similarity matrix

   python MakeSimMatrix.py myhomopolymers 4 mySimMatrix.json

   ** Make a note of the number of neighbours listed in the output **

   