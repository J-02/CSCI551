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
    # getting
    C, OCC, bwtSeq, seq = fmindex(file, Print=False)
    top = 0
    bottom = len(bwtSeq) - 1

    for char in reversed(Q):
        if char in C:
            top = C[char] + (OCC[char, top-1] if top > 0 else 0)
            bottom = C[char] + OCC[char, bottom] - 1
        else:
            # Character not in text, return empty range
            return None, None

    print(f'range[{seq}, {Q}] = [{top}, {bottom}]')

    return top, bottom + 1

def suffixArray(string):
    suffixes = []
    for idx, i in enumerate(string):
        suffixes.append(string[idx:])
    SA = sorted(suffixes)
    return SA

def sampledSuffixArray(s, k=10):
    SA = suffixArray(s)
    sampledSA = {i: pos for i, pos in enumerate(SA) if i % k == 0}
    return sampledSA, k


def findEntry(seq, bwt_seq, C, OCC, sampled_suffix_array, k, target_index):
    if target_index in sampled_suffix_array:
        return sampled_suffix_array[target_index]

    # Find the closest sampled index below the target index
    sampled_index = (target_index // k) * k

    # Initialize the current index and character to the sampled index and corresponding character
    current_index = sampled_index
    current_char = seq[sampled_index]

    # Follow the FM-index backward until reaching the target index
    while current_index != target_index:
        # Compute the previous index using the FM-index
        previous_index = C[current_char] + OCC[current_char, current_index - 1]

        # Update the current index and character
        current_index = previous_index
        current_char = bwt_seq[current_index]

    # Return the position in seq corresponding to the target index
    return sampled_suffix_array[sampled_index] + (target_index - sampled_index)



print("Part 1:\n S = ACTGGGAAATCGAAGACCCGG")
file = "test.fa"
fmindex(file)

print("\nPart 2:\n")

test1 = getSeq(file)
print(f"Suffix array for: {test1}")
print(suffixArray(test1))
range_search(file, Q="AT")
range_search(file, Q="AA")
file = "WedTest.fa"
test2 = getSeq(file)
print(f"\nSuffix array for: {test2}")
print(suffixArray(test1))
range_search(file, Q="AT")
range_search(file, Q="AA")
file = "test2.fa"
test3 = getSeq(file)
print(f"\nSuffix array for: {test3}")
print(suffixArray(test1))
range_search(file, Q="AT")
range_search(file, Q="AA")

print("\nPart 3:\n")
file = "test3.fa"
test3 = getSeq(file)

sSA, k = sampledSuffixArray(test3)
C, OCC, bwtSeq, seq = fmindex(file, Print=False)
entry = findEntry(seq, bwtSeq, C, OCC, sSA, k, 5)
print(entry)

