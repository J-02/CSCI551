import sys
import numpy as np
# runs in terminal with "python main.py test.fasta search"
# python main.py test.fasta TAC
# can be run on any .fasta file in directory
# This code reads the fasta file to dictionary

def main(scores, file):
    # Search each sequence in file for string and return z array and exact pattern matching
    motifs = readFasta(file)

    T = "-" + motifs[1]
    S = "-" + motifs[0]
    print(T)
    print(S)
    V = np.zeros(shape = (len(T), len(S)))
    print(V)



def readFasta(file):
    motifs = []
    with open(file) as f:
        lines = f.readlines()
    for i in range(0, len(lines)):
        seq = lines[i].strip()  # skips the line if > is present
        if seq[0] == '>':  # Executes check of '>' character
            text = ""
        else:
            motifs.append(seq)

    print(motifs)
    return motifs

main((2,-1,-1), "test.fa")

# to run in terminal
#if __name__ == "__main__":
#   a = int(sys.argv[1])
#   b = int(sys.argv[2])
#   c = int(sys.argv[3])
#   file = sys.argv[4]
#   scoring = (a,b,c)
#   main(scoring, file)
