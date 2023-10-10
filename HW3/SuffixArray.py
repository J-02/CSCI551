import sys
import numpy as np
# runs in terminal with "python SuffixArray.py test.fa CT"
# python SuffixArray.py test.fa CT
# search term optional, returns suffix array otherwise
# can be run on any text file in directory with 1 sequence

def main(file, Q=None):
    motif = readFasta(file)[0]
    SA = suffixArray(motif)
    for idx, i in enumerate(SA):
        print(f'SA[{idx}] = {i}')
    if Q:
        return binarySearch(Q,SA)
    return SA

def suffixArray(string):
    suffixes = []
    for idx, i in enumerate(string):
        suffixes.append(string[idx:])
    SA = sorted(suffixes)
    return SA

def binarySearch(Q, SA):
    # setting initial L and R
    L = 0
    R = len(SA)
    # looping through search
    while True:
        # getting middle index M
        M = (L + R) // 2
        # getting longest common prefix m
        m = lcp(Q,SA[M])

        # checks if length of the longest common prefix equals the length of the pattern
        if len(m) == len(Q):
            print(f"Match found: SA[{M}] = {m}")
            return SA[M]

        #if middle string lexographically after search term
        elif SA[M][:-1] > Q:
            print(SA[M][:-1], Q)
            R = M
            continue
        # SA[m] occurs before Search pattern
        else:
            # if we have gone through all of SA
            if L == R:
                print("Match not found")
                return None
            # changes pointer
            else:
                if (L + 1) != R:
                    L = M
                else:
                    L = R

def lcp(s1, s2):
    len1, len2 = len(s1), len(s2)
    prefix = []

    i, j = 0, 0
    while i < len1 and j < len2:
        if s1[i] == s2[j]:
            prefix.append(s1[i])
            i += 1
            j += 1
        else:
            break

    # Convert the list of common characters into a string and store it in 's'
    s = ''.join(prefix)
    return s


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
# main("test.fa", search)

# to run in terminal
if __name__ == "__main__":
    file = sys.argv[1]
    Q = sys.argv[2]
    main(file, Q)

