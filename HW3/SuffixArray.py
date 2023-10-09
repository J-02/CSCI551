import sys
import numpy as np
# runs in terminal with "python myAlign.py 2 -1 -1 test.fa"
# python myAlign.py 2 -1 -1 test.fa
# can be run on any text file in directory with 2 sequences
# This code reads the text file to list of 2 elements

def main(scores, file):

    motif = readFasta(file)

def suffixArray(string):
    return


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
            motifs.append(seq)

    print(f"Getting Global alignment for {motifs}")
    return motifs

# test run function
# scores in tuple (a, b, c)
# main((2,-1,-1), "test.fa")

# to run in terminal
if __name__ == "__main__":
   a = int(sys.argv[1])
   b = int(sys.argv[2])
   c = int(sys.argv[3])
   file = sys.argv[4]
   scoring = (a,b,c)
   main(scoring, file)
