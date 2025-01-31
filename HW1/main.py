import sys
#exact pattern match utilizing the Z algorithm
# runs in terminal with "python main.py test.fasta search"
# python main.py test.fasta TAC
# can be run on any .fasta file in directory (with single and multiple sequences)
# This code reads the fasta file to dictionary = 'motif'

def main(file, search):
    # Search each sequence in file for string and return z array and exact pattern matching
    motifs = readFasta(file)

    for key, seq in motifs.items():
        s = search+"#"+seq  # creates concatenated string with search term
        zarray = zalg(seq)  # zarray without search term
        exact = zalg(s)  # zarray with search term
        indices = [i for i, v in enumerate(exact) if v == len(search)]  # indices where exact pattern starts
        print(f"Searching {key}:{seq} for {search}")
        print(f"Z-array: {zarray}")
        print(f"Exact pattern matching starting at indicies: {indices}\n\n")

def readFasta(file):
    motifs = dict()
    key = ""
    with open(file) as f:
        lines = f.readlines()
    for i in range(0, len(lines)):
        seq = lines[i].strip()  # skips the line if > is present
        if seq[0] == '>':  # Executes check of '>' character
            if key:
                motifs[key] = text
            key = seq[
                  1:]  # Following lines are store into the key, which is where the sequences in the fasta will be stored
            text = ""
        else:
            text += seq
    motifs[key] = text
    return motifs


def zalg(s):
    # returns Z array length of input sequence
    n = len(s)
    Z = [0] * n  # empty z-array the length of sequence
    Z[0] = n
    l = r = 0

    for k in range(1, n):
#iterate through all k values in range n-1
        if k <= r: #Case1: k is within the interval [1,r]
            Z[k] = min(r - k + 1, Z[k - l]) #determkne the distance from right most boundatry from r to k and find the Z value at the position k-1
#Case 2: execute explicit charactet comparison to calculate Z[k]
        while k + Z[k] < n and s[Z[k]] == s[k + Z[k]]: # within the bounds of s and characters match, increment Z[k]+1 to the length of the matching sequence
            Z[k] += 1
#Update the interval to the extent of the curren Z[k]
        if k + Z[k] - 1 > r:
            l, r = k, k + Z[k] - 1
#Return calculated Z value
    return Z

# to run in terminal
if __name__ == "__main__":
    file = (sys.argv[1])
    search = (sys.argv[2])
    main(file, search)
