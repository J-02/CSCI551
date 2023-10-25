from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
def getSeq(file):
    record = SeqIO.read(file, "fasta")
    seq = Seq(str(record.seq)+ "$")
    return seq

def bwt(s):
    # make sure s is a string
    s = str(s)

    # create a list/table of all cyclic rotations of 's' then sort them lexographically
    bw = sorted(s[i:] + s[:i] for i in range(len(s)))

    # take the last column to get the BWT
    last_column = [row[-1] for row in bw]

    # return the string
    return "".join(last_column)

def fmindex(file, Print=True):
    # reading the file into string format
    seq = getSeq(file)

    # getting the bwt of the seq
    bwtSeq = bwt(seq)

    # Compute the C table as a dictionary
    C = {}
    # get occurrences of each unique letter in sequence
    charC = Counter(bwtSeq)
    # initializing charSum
    charSum = 0
    # sort char counts lexicographically then
    for char in sorted(charC.keys()):
        C[char] = charSum
        charSum += charC[char]
        # gets the amount of characters lexicographically smaller than char

    # Compute the occurrence table (OCC)
    OCC = {}
    # initializing empty char counts
    char_counts = {char: 0 for char in set(bwtSeq)}

    # loop through all the chars and indexes in the BWT Seq
    for i, char in enumerate(bwtSeq):
        # count the occurrence
        char_counts[char] += 1

        # loop through table and update where the occurrences occur
        for char, count in char_counts.items():
            # update count
            OCC[char, i] = count
    if not Print:
        return C, OCC, bwtSeq, seq
    # Output the FM-index
    print(f'BW = {bwtSeq}')
    for char, count in C.items():
        print(f'C[{char}] = {count}')
    for (char, i), count in OCC.items():
        print(f'OCC[{char},{i}] = {count}')

    return C, OCC, bwtSeq, seq


def range_search(file, Q):
    C, OCC, bwtSeq, seq =fmindex(file, Print=False)
    top = 0
    bottom = len(bwtSeq) - 1

    for char in reversed(Q):
        if char in C:
            top = C[char] + (OCC[char, top - 1] if top > 0 else 0)
            bottom = C[char] + OCC[char, bottom] - 1
        else:
            # Character not in text, return empty range
            return None, None

    print(f'range[{seq}, {Q}] = [{top}, {bottom+1}]')

    return top, bottom + 1

fmindex("test.fa")
range_search("test.fa", Q="AT")
