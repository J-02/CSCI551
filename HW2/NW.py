import sys
import numpy as np

# Function to read sequences from a FASTA file (already provided)
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
    return motifs

def main(scores, file):
    # Read sequences from the FASTA file using your existing function
    motifs = readFasta(file)

    # Extract sequences and name them accordingly (S and T)
    S = "-" + motifs[0]  # Adding a gap to the beginning
    T = "-" + motifs[1]  # Adding a gap to the beginning

    # Extract scoring parameters
    match = scores[0]
    mismatch = scores[1]
    insertion = scores[2]
    deletion = scores[3]
    gap_penalty = -1  # Gap penalty (you can customize this)

    # Initialize the scoring matrix V
    m = len(T)
    n = len(S)
    V = np.zeros(shape=(m + 1, n + 1))

    # Initialize the first row and first column of V with gap penalties
    for i in range(1, m + 1):
        V[i][0] = i * gap_penalty
    for j in range(1, n + 1):
        V[0][j] = j * gap_penalty
###
    # Filling in the scoring matrix V
    for i in range(1, m):
        for j in range(1, n):
            match_score = V[i - 1][j - 1] + (match if T[i-1] == S[j-1] else mismatch)
            gap_in_T = V[i - 1][j] + insertion
            gap_in_S = V[i][j - 1] + deletion
            V[i][j] = max(match_score, gap_in_T, gap_in_S)

###
    # Filling in the scoring matrix V
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_score = V[i - 1][j - 1] + (match if T[i - 1] == S[j - 1] else mismatch)
            gap_in_T = V[i - 1][j] + insertion
            gap_in_S = V[i][j - 1] + deletion
            V[i][j] = max(match_score, gap_in_T, gap_in_S)

    # The optimal global alignment score is in V[m][n]
    optimal_score = V[m][n]

    # Print the optimal alignment score (you can add alignment retrieval here)
    print("Optimal Alignment Score:", optimal_score)

# Test run function with example scores and file
main((2, -6, -6, -6), "test.fa")

# To run in the terminal, uncomment and modify the following lines:
# if __name__ == "__main__":
#     a = int(sys.argv[1])
#     b = int(sys.argv[2])
#     c = int(sys.argv[3])
#     d = int(sys.argv[4])
#     file = sys.argv[5]
#     scoring = (a, b, c, d)
#  main(scoring, file)
