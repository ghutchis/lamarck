#!/bin/python

# The goal of this script is to randomly choose n lines from a text files by line without duplicate

import sys, os, pdb, random

#The prefix will add to the front of each line
prefix = ""

if len(sys.argv) != 3 and len(sys.argv) != 4:
    print("syntax: randomlines.py <input_file> <num_of_sample>")
    sys.exit()
else:
    input_file = sys.argv[1]
    num_of_sample = int(sys.argv[2])

if len(sys.argv) == 4:
    prefix = sys.argv[3]

#pdb.set_trace()
lines = open(input_file, "r").read().splitlines()
if len(lines) < num_of_sample:
    print("The number of random lines you asked for is larger than the lines in the target file.\n" +
        "Shrink it to lines of target file automatically.\n" +
        "-----------------------------------------------\n")
    num_of_sample = len(lines)

chosen_lines_num = random.sample(range(0,len(lines)), num_of_sample)   #range(0,3) only generate 0,1,2
for i in chosen_lines_num:
    print(prefix + lines[i])