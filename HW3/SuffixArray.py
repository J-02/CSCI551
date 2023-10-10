import sys
import numpy as np
# runs in terminal with "python myAlign.py 2 -1 -1 test.fa"
# python myAlign.py 2 -1 -1 test.fa
# can be run on any text file in directory with 2 sequences
# This code reads the text file to list of 2 elements

def main(file):
    motif = readFasta(file)[0]
    SA = suffixArray(motif)
    for idx, i in enumerate(SA):
        print(f'SA[{idx}] = {i}')

def suffixArray(string):
    suffixes = []
    for idx, i in enumerate(string):
        suffixes.append(string[idx:])
    SA = sorted(suffixes)
    return SA

def readFasta(file):
    motifs = []
    with open(file) as f:
        lines = f.readlines()
    for i in range(0, len(lines)):
        seq = lines[i].strip()  # skips the line if > is present
        if seq == "": # if there is a space
            continue
        elif seq[0] == '>':  # Executes check of '>' character
            text = ""

        else:
            motifs.append(seq+"$")

    print(f"Getting Suffix Array for {motifs}")
    return motifs

# test run function
# scores in tuple (a, b, c)
# main((2,-1,-1), "test.fa")

# to run in terminal
if __name__ == "__main__":
   file = "test.fa"
   main(file)
