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

def fmindex(file):
    # reading the file into string format
    seq = getSeq(file)

    # getting the bwt of the seq
    bwt_seq = bwt(seq)

    # Compute the C table
    C = {}
    char_counts = Counter(bwt_seq)
    char_sum = 0
    for char in sorted(char_counts.keys()):
        C[char] = char_sum
        char_sum += char_counts[char]

    # Compute the occurrence table (OCC)
    OCC = {}
    char_counts = {char: 0 for char in set(bwt_seq)}
    for i, char in enumerate(bwt_seq):
        char_counts[char] += 1
        for char, count in char_counts.items():
            OCC[char, i] = count

    # Output the FM-index
    print(f'BW = {bwt_seq}')
    for char, count in C.items():
        print(f'C[{char}] = {count}')
    for (char, i), count in OCC.items():
        print(f'OCC[{char},{i}] = {count}')

fmindex("test.fa")
