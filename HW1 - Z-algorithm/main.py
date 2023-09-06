import sys

def main(file, search):
    # Search each sequence in file for string and return z array and exact pattern matching
    motifs = readFasta(file)

    for key, seq in motifs.items():
        s = search+"#"+seq
        zarray = zalg(seq)
        exact = zalg(s)[len(search)+1:]
        indices = [i for i, v in enumerate(exact) if v == len(search)]
        print(f"Searching {key}:{seq} for {search}")
        print(f"Z-array: {zarray}")
        print(f"Exact pattern matching at indicies: {indices}")

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
    ## returns Z array length of input sequence
    n = len(s)
    Z = [0] * n
    Z[0] = n

    l = r = 0
    for i in range(1, n):
        if i <= r:
            Z[i] = min(r - i + 1, Z[i - l])

        while i + Z[i] < n and s[Z[i]] == s[i + Z[i]]:
            Z[i] += 1

        if i + Z[i] - 1 > r:
            l, r = i, i + Z[i] - 1

    return Z

if __name__ == "__main__":
    file = (sys.argv[1])
    search = (sys.argv[2])
    main(file, search)
