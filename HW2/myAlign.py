import sys
import numpy as np
# runs in terminal with "python myAlign.py 2 -1 -1 test.fa"
# python myAlign.py 2 -1 -1 test.fa
# can be run on any text file in directory with 2 sequences
# This code reads the text file to list of 2 elements

def main(scores, file):

    motifs = readFasta(file)
    # get S and T
    match = scores[0]
    # score for match
    mismatch = scores[1]
    # score for mismatch
    insdel = scores[2]
    # score for insertion/deletion

    gap_penalty = insdel


    # adding gaps to the beginning of each sequence
    S = "-" + motifs[1]
    T = "-" + motifs[0]

    m = len(T)-1
    n = len(S)-1
    # initialize V
    V = np.zeros(shape=(m+1, n+1))

    # filling first row and column with gap penalties
    for i in range(n+1):
        V[0,i] = i*gap_penalty
    for i in range(m+1):
        V[i,0] = i*gap_penalty

    # iterating through each letter in T not including gap
    for i in range(1,m+1):
        # i is index of Ti in T
        # iterating through each letter in S not including gap
        for j in range(1,n+1):
            # j is index of Sj in S
            if T[i] == S[j]:
                # if equal add match score
                V[i,j] = match + V[i-1,j-1]
            elif T[i] != S[j]:
                # checking others if not a match
                V[i,j] = max(V[i][j-1] + gap_penalty, V[i-1][j] + gap_penalty, V[i-1][j-1] + mismatch)
    print(f"Optimal Global Alignment: {int(V[m,n])}")

    all_alignments = []
    # print(V)
    # recursive_opt_alignment recursively calls itself to get all optimal global alignments from V by tracing the values in it
    def recursive_opt_alignment(i, j, aT, aS, alignments):
        # adds full optimal alignments to alignments
        if i == 0 and j == 0:
            alignments.append((''.join(reversed(aT)), ''.join(reversed(aS))))
            return

        # continues the global alignment from the diagonal
        if i > 0 and j > 0:
            diagonal = V[i - 1][j - 1]
            if V[i][j] == diagonal + (match if T[i] == S[j] else mismatch):
                recursive_opt_alignment(i - 1, j - 1, aT + [T[i]], aS + [S[j]], alignments)

        # continues the global alignment from above V[i,j] at V[i-1,j]
        if i > 0:
            vertical = V[i - 1][j]
            if V[i][j] == vertical + gap_penalty:
                recursive_opt_alignment(i - 1, j, aT + [T[i]], aS + ['-'], alignments)

        # continues the global alignment from next to V[i,j] to V[i,j-1]
        if j > 0:
            horizontal = V[i][j - 1]
            if V[i][j] == horizontal + gap_penalty:
                recursive_opt_alignment(i, j - 1, aT + ['-'], aS + [S[j]], alignments)

    recursive_opt_alignment(m, n, [], [], all_alignments)

    print(f"All optimal global alignments {all_alignments}")
    # Reverse the aligned sequences to get the final alignment


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
