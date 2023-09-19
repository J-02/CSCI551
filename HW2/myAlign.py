import sys

# runs in terminal with "python main.py test.fasta search"
# python main.py test.fasta TAC
# can be run on any .fasta file in directory
# This code reads the fasta file to dictionary

def main(file, search):
    # Search each sequence in file for string and return z array and exact pattern matching
    motifs = readFasta(file)

    T = motifs[1]
    S = motifs[0]



def readFasta(file):
    motifs = []
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
    motifs.append(text)
    return motifs


# to run in terminal
if __name__ == "__main__":
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    c = int(sys.argv[3])
    file = sys.argv[4]
    scoring = (a,b,c)
    main(scoring, file)
