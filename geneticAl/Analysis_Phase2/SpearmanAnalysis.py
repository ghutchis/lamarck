#!/usr/bin/python
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from optparse import OptionParser


# Read in output file from SpearmanCorrelations.py script.
dataFile = sys.argv[1]
df = pd.read_csv(dataFile, sep='\t', header=None)
df.columns =  ['monomers','generations','percentile','sp_corr']

def analysis(perc):
    # Sort data to use given percentile
    perc_sorted = df[(df.percentile == perc)]

    print perc_sorted

    # Sort each group of monomers into its own variable
    perc_131 = perc_sorted[(perc_sorted.monomers == 131)]
    perc_442 = perc_sorted[(perc_sorted.monomers == 442)]
    perc_611 = perc_sorted[(perc_sorted.monomers == 611)]
    perc_909 = perc_sorted[(perc_sorted.monomers == 909)]
    perc_1235 = perc_sorted[(perc_sorted.monomers == 1235)]
    perc_1759 = perc_sorted[(perc_sorted.monomers == 1759)]

    # Set x and y for each data set for ease in plotting
    x1 = perc_131.generations
    y1 = perc_131.sp_corr
    x2 = perc_442.generations
    y2 = perc_442.sp_corr
    x3 = perc_611.generations
    y3 = perc_611.sp_corr
    x4 = perc_909.generations
    y4 = perc_909.sp_corr
    x5 = perc_1235.generations
    y5 = perc_1235.sp_corr
    x6 = perc_1759.generations
    y6 = perc_1759.sp_corr

    # Generate subplot with each data set as its own scatter plot
    f, axarr = plt.subplots(3, 2)
    f.set_tight_layout(True)

    #f.tight_layout(pad=0.5, w_pad=0.5, h_pad=1.0)
    axarr[0, 0].scatter(x1, y1, s=30, alpha=0.15, marker='o')
    axarr[0, 0].set_title('131 Monomers')
    axarr[0, 0].set_xlabel('Generations')
    axarr[0, 0].set_ylabel('Spearman Correlation')
    axarr[0, 1].scatter(x2, y2, s=30, alpha=0.15, marker='o')
    axarr[0, 1].set_title('442 Monomers')
    axarr[1, 0].scatter(x3, y3, s=30, alpha=0.15, marker='o')
    axarr[1, 0].set_title('611 Monomers')
    axarr[1, 1].scatter(x4, y4, s=30, alpha=0.15, marker='o')
    axarr[1, 1].set_title('909 Monomers')
    axarr[2, 0].scatter(x5, y5, s=30, alpha=0.15, marker='o')
    axarr[2, 0].set_title('1235 Monomers')
    # Set [2,1] as blank plot because have only 5 data sets
    axarr[2,1].scatter(x6, y6, s=30, alpha=0.15, marker='o')
    axarr[2, 0].set_title('1759 Monomers')
    # Save plot
    savefig('ScatterPlots-%s.pdf' % perc)
    # Show plot gives error...don't know why
    # plt.show()

    perc_sorted.boxplot(column='sp_corr', by=['monomers', 'generations'], notch=False)
    # set your own proper title
    plt.title("Box Plots of All Data")
    # get rid of the automatic 'Boxplot grouped by group_by_column_name' title
    plt.suptitle("")
    savefig('BoxPlots-%s.pdf' % perc, bbox_inches='tight')

    perc_131.boxplot(column='sp_corr', by=['monomers', 'generations'], notch=False)
    # set your own proper title
    plt.title("131 Monomers-%s" % perc)
    # get rid of the automatic 'Boxplot grouped by group_by_column_name' title
    plt.suptitle("")
    savefig('BoxPlots-131-%s.pdf' % perc, bbox_inches='tight')
    perc_442.boxplot(column='sp_corr', by=['monomers', 'generations'], notch=False)
    plt.title("442 Monomers-%s" % perc)
    savefig('BoxPlots-442-%s.pdf' % perc)
    perc_611.boxplot(column='sp_corr', by=['monomers', 'generations'], notch=False)
    plt.title("611 Monomers-%s" % perc)
    savefig('BoxPlots-611-%s.pdf' % perc)
    perc_909.boxplot(column='sp_corr', by=['monomers', 'generations'], notch=False)
    plt.title("909 Monomers-%s" % perc)
    savefig('BoxPlots-909-%s.pdf' % perc)
    perc_1235.boxplot(column='sp_corr', by=['monomers', 'generations'], notch=False)
    plt.title("1235 Monomers-%s" % perc)
    savefig('BoxPlots-1235-%s.pdf' % perc)
    perc_1759.boxplot(column='sp_corr', by=['monomers', 'generations'], notch=False)
    plt.title("1759 Monomers-%s" % perc)
    savefig('BoxPlots-1759-%s.pdf' % perc)

    # #fit = np.polyfit(a+b*log(c*x),y,1)
    # #fit_fn = poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    #
    # par = np.polyfit(x1, y1, 1, full=True)
    # slope=par[0][0]
    # intercept=par[0][1]
    # xl = [min(x1), max(x1)]
    # yl = [slope*xx + intercept  for xx in xl]
    #
    # # coefficient of determination, plot text
    # variance = np.var(y1)
    # residuals = np.var([(slope*xx + intercept - yy)  for xx,yy in zip(x1,y1)])
    # Rsqr = np.round(1-residuals/variance, decimals=2)
    # plt.text(.9*max(x1)+.1*min(x1),.9*max(y1)+.1*min(y1),'$R^2 = %0.2f$'% Rsqr, fontsize=30)
    #
    # #plt.xlabel("X Description")
    # #plt.ylabel("Y Description")
    #
    # # error bounds
    # yerr = [abs(slope*xx + intercept - yy)  for xx,yy in zip(x1,y1)]
    # parerr = np.polyfit(y1, yerr, 2, full=True)
    #
    # yerrUpper = [(xx*slope+intercept)+(parerr[0][0]*xx**2 + parerr[0][1]*xx + parerr[0][2]) for xx,yy in zip(x1,y1)]
    # yerrLower = [(xx*slope+intercept)-(parerr[0][0]*xx**2 + parerr[0][1]*xx + parerr[0][2]) for xx,yy in zip(x1,y1)]
    #
    # plt.plot(xl, yl, '-r')
    # plt.plot(x1, yerrLower, '--r')
    # plt.plot(x1, yerrUpper, '--r')
    # plt.show()

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
