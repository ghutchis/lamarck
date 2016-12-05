#!/usr/bin/python
import sys
import pandas as pd
import numpy as np


def analysis(perc):
    # Makes a pivot table with SMILES as the rows and columns of different row numbers
    pt = pd.pivot_table(df, index=["SMILES"], values=["Run"], columns=["Run"], fill_value=0)
    # Add across each row to get the total number of monomers and sort from largest to smallest
    pt_sum = pt.sum(1).order(ascending=False)

    # Takes the values that are larger than the quartile being analyzed
    a = pt_sum[pt_sum > pt_sum.quantile(perc)]

    # Identifies the number of monomers within that top percentage which will be used in the analysis
    num_monomers = len(a.index)
    print num_monomers

    # Generates smaller pivot table which includes only the top monomers
    pt_subset = pt.loc[a.index.values]

    # Pairwise spearman correlation of all runs
    corr = pt_subset.corr(method='spearman')
    # Average spearman correlation ontained by removing duplicate calculations
    pt_avg = corr.values[np.triu_indices_from(corr.values, 1)].mean()

    # Calculate std error by making duplicate averages nan and the taking stdev/ sqrt(number of calcs)
    c = corr.copy()
    c.values[np.tril_indices_from(c)] = np.nan
    s = c.unstack().mean()
    s_stderr = corr.values[np.triu_indices_from(corr.values, 1)].std() / np.sqrt(s.count())
    # Saves info to output file 
    #output_file.write("Percentage=\t%s\tNumber of monomers =\t %s\t Corr=\t %.3f\tStdErr =\t%.3f\n %s\n"
     #                 % (perc, num_monomers, pt_avg, s_stderr,a.index.values))
    output_file.write("Percentage=\t%s\tNumber of monomers =\t %s\t Corr=\t %.3f\tStdErr =\t%.3f\n"
                      % (perc, num_monomers, pt_avg, s_stderr))

# Read in output file from MonomerCounts.py script. Make sure there are headers on the column named Run, SMILES, Count
dataFile = sys.argv[1]
df = pd.read_csv(dataFile, sep='\t', header=0)
df.columns = ['Run', 'SMILES', 'Counts']
output_file = open(sys.argv[2], "w")

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

output_file.close()
