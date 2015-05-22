#!/usr/bin/python
import sys
import pandas as pd
import numpy as np

# To run this file:
# Argument 1: Output file from MonomerCounts.py script with column headers Run, SMILES, Count.
# Argument 2: File name to which the program should save the output. Script will append to existing file.
# Argument 3: Number of generations in the file.
# Argument 4: Number of monomers in the data set.


def analysis(perc):
    # Makes a pivot table with SMILES as the rows and columns of different row numbers
    pt = pd.pivot_table(df, index=["SMILES"], values=["Run"], columns=["Run"], fill_value=0)
    # Add across each row to get the total number of monomers and sort from largest to smallest
    pt_sum = pt.sum(1).order(ascending=False)

    # Takes the values that are larger than the quartile being analyzed
    a = pt_sum[pt_sum > pt_sum.quantile(perc)]

    # Identifies the number of monomers within that top percentage which will be used in the analysis
    num_monomers = len(a.index)

    # Generates smaller pivot table which includes only the top monomers
    pt_subset = pt.loc[a.index.values]

    # Pairwise spearman correlation of all runs
    corr = pt_subset.corr(method='spearman')

    # Removes duplicate Spearman correlations and substituting with nan
    c = corr.copy()
    c.values[np.tril_indices_from(c)] = np.nan

    s = c.unstack().mean()

    # Generate ndarray of values of the Spearman correlations
    s1 = s.values

    # Removes nan values leaving just the desired  Spearman correlations
    sv = s1[~np.isnan(s1)]

    # Converts numpy array to list of values which can then be saved individually
    sv_values = list(sv)

    # Saves the values in a test file in the format percentage, Generations, Spearman Value which
    # can then be plotted. This file will just be appended to which allows many data sets to be combined into
    # one long list with all data which can then be plotted, analyzed, etc.

    for i in sv_values:
        generations = int(sys.argv[3])
        monomers = int(sys.argv[4])
        with open(sys.argv[2], "a") as output_file:
            output_file.write("%i\t%i\t%s\t%f\n" % (monomers, generations, perc, i))



# Read in output file from MonomerCounts.py script. Make sure there are headers on the column named Run, SMILES, Count
dataFile = sys.argv[1]
df = pd.read_csv(dataFile, sep='\t', header=0)



# What percentiles to we want to return? When 0.75 is run, the output is the top 25% of the data, etc.

perc1 = 0.75
analysis(perc1)

perc2 = 0.80
analysis(perc2)

perc3 = 0.85
analysis(perc3)

perc4 = 0.90
analysis(perc4)

perc5 = 0.95
analysis(perc5)


