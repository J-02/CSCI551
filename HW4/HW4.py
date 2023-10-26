from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

def getSeq(file):
    record = SeqIO.read(file, "fasta")
    seq = Seq(str(record.seq)+ "$")
    return seq

def bwt(s):
    # create a list/table of all cyclic rotations of 's' then sort them lexographically
    bw = sorted(s[i:] + s[:i] for i in range(len(s)))
    # take the last column to get the BWT
    last_column = [row[-1] for row in bw]
    # return the string concatenated string
    return "".join(last_column)

def suffixArray(s):
    # takes each suffix and place and puts suffix index in lex sorted list
    return [i for _, i in sorted((s[i:], i) for i in range(len(s)))]


class Sequence:
    def __init__(self, file, k=10):
        # initializes
        # k: for the sampled suffix array
        # seq: to string with '$' at the end
        # bwtSeq: the BWT sequence
        # C, OCC: C and OCC arrays from fmindex
        # SA: suffix array to check correctness of sampled SA and to construct SSA from
        # SSA: sampled suffix array with 'k' sampling rate
        # suffixes: list of suffixes to find by index
        self.k = k
        self.seq = getSeq(file)
        self.bwtSeq = bwt(self.seq)
        self.C, self.OCC = self.fmindex()
        self.SA = suffixArray(self.seq)
        self.SSA = self.sampledSuffixArray(k=self.k)
        self.suffixes = self.getSuffixes()
        self.len = len(self.bwtSeq)

    def getSuffixes(self):
        # gets dictionary of index and suffix not sorted
        suffixes = {(idx, self.seq[idx:]) for idx, i in enumerate(self.seq)}
        return suffixes

    def printFmindex(self):
        # prints FM index information
        print(f"Sequence: {self.seq}")
        print(f'BW = {self.bwtSeq}')
        for char, count in self.C.items():
            print(f'C[{char}] = {count}')
        for (char, i), count in self.OCC.items():
            print(f'OCC[{char},{i}] = {count}')
        return

    def fmindex(self):
        # gets the C and OCC arrays
        C = {}
        # get occurrences of each unique letter in sequence
        charC = Counter(self.bwtSeq)
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
        char_counts = {char: 0 for char in set(self.bwtSeq)}

        # loop through all the chars and indexes in the BWT Seq
        for i, char in enumerate(self.bwtSeq):
            # count the occurrence
            char_counts[char] += 1

            # loop through table and update where the occurrences occur
            for char, count in char_counts.items():
                # update count
                OCC[char, i] = count
        return C, OCC

    def printRangeSearch(self, Q):
        # prints the range search information
        top, bottom = self.rangeSearch(Q)
        print(f'range[{self.seq}, {Q}] = [{top}, {bottom-1}]')

    def rangeSearch(self, Q):
        # performs range search for string Q
        top = 0
        bottom = len(self.bwtSeq) - 1

        for char in reversed(Q):
            if char in self.C:
                top = self.C[char] + (self.OCC[char, top - 1] if top > 0 else 0)
                bottom = self.C[char] + self.OCC[char, bottom] - 1
            else:
                # Character not in text, return empty range
                return None, None

        return top, bottom + 1

    def sampledSuffixArray(self, k=10):
        # creates samples suffix array from SA with k rate
        sampledSA = {i: pos for i, pos in enumerate(self.SA) if i % k == 0}
        return sampledSA

    def printFindEntry(self,i):
        # prints the info for finding entry
        print(f"SA[{i}] = {self.findEntry(i)} ")

    def findEntry(self, i):
        i_prime = i
        delta = 0
        while i_prime not in self.SSA:
            x = self.bwtSeq[i_prime]
            i_prime = self.C[x] + self.OCC[(x, i_prime)] - 1
            delta += 1
        return (self.SSA[i_prime] + delta) % self.len


seq1 = Sequence('test.fa')
seq2 = Sequence('test2.fa')
seq3 = Sequence('WedTest.fa')
seq4 = Sequence('test3.fa', 250)


print("Part 1:\nFM-Index:")
seq1.printFmindex()
print()

query = 'AT'
print(f"Part 2:\n Range search:")
seq1.printRangeSearch('AT')
seq1.printRangeSearch('AA')
seq2.printRangeSearch('AT')
seq2.printRangeSearch('AA')
seq3.printRangeSearch('AT')
seq3.printRangeSearch('AA')
print()

# Answer to Bonus question: Suffix Array Sampling
print("Part 3:\n Finding SA[i] using sampled suffix array where k = 250 for 'test3.fa':")
print(f"Length of test4.fa: {seq4.len}")
print(f"Samples Suffix Array:\n {seq4.SSA}")
i = 34
print(f"Finding entry {i}:")
seq4.printFindEntry(i)
print(f"From actual SA[{i}]: {seq4.SA[i]}")
