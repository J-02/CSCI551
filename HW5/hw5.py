from Bio import SeqIO
from Bio import Align
def read_fasta(file):
    # reading the fasta seqs
    sequences = {}
    for record in SeqIO.parse(file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def center_star(seqs, alpha, beta, output_center_sequence=False):

    def delta(x, y):
        # sets the scoring
        return 0 if x == y else alpha if y != '-' else beta

    def seq_distance(seq1, seq2):
        # gets distance between a pair
        distance = 0
        len_seq1, len_seq2 = len(seq1), len(seq2)
        for i in range(max(len_seq1, len_seq2)):
            char_seq1 = seq1[i] if i < len_seq1 else '-'
            char_seq2 = seq2[i] if i < len_seq2 else '-'
            distance += delta(char_seq1, char_seq2)
        return distance
    def pairwise_distances():
        # totals up distances for each sequence
        distances = {}
        names = list(seqs.keys())
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                seq1, seq2 = names[i], names[j]
                distances[(seq1,seq2)] = seq_distance(seqs[seq1],seqs[seq2])
        return distances

    def center_sequence(distances):
        # finds the center star by getting min distance
        total_distances = {seq: 0 for seq in seqs}
        for (seq1, seq2), distance in distances.items():
            total_distances[seq1] += distance
            total_distances[seq2] += distance
        return min(total_distances, key=total_distances.get), total_distances

    distances = pairwise_distances()
    center, total_distances = center_sequence(distances)

    if output_center_sequence:
        for i in seqs.keys():
            print(f"{i}: {seqs[i]}")

        print(f"Scoring: Alpha = {alpha}, Beta = {beta}")
        print("Center Sequence:", center, seqs[center])
        print("Total Distances:", total_distances)
    else:
        return center

def align_to_center(center, sequence):
    # Aligns given seq to the center
    alignments = aligner.align(center, sequence)
    best_alignment = alignments[0]
    # Extracting aligned sequences
    test = best_alignment.format()
    lines = test.splitlines()
    center = lines[0].split()[2]  # Extracting the sequence part
    sequence = lines[2].split()[2]
    return center, sequence

def find(s, ch):
    # finds where the new gaps are
    return [i for i, ltr in enumerate(s) if ltr == ch]

def find_new_gaps(old_center, new_center):
    # Find the new gaps
    old_gaps = set(find(old_center, "-"))
    new_gaps = set(find(new_center, "-"))
    return list(new_gaps - old_gaps)

def MSA(center_seq, seqs):
    # getting the MSA by aligning each seq with center and entering gaps when the seq is larger than the center
    center = center_seq
    aligned_seqs = {}
    for seq_key in seqs.keys():
        center_seq2, sequence = align_to_center(center, seqs[seq_key])
        if center_seq2 != center:
            gaps = find_new_gaps(center, center_seq2)
            for gap_index in gaps:
                for aligned_key in aligned_seqs.keys():
                    aligned_seqs[aligned_key] = aligned_seqs[aligned_key][:gap_index] + "-" + aligned_seqs[aligned_key][gap_index:]
            center = center_seq2
        aligned_seqs[seq_key] = sequence

    return center, aligned_seqs



alpha = 1
beta = 1
fasta = "test.fasta"  # Replace with your file path
sequences = read_fasta(fasta)
center_star(sequences, alpha, beta, True)
center = center_star(sequences, alpha, beta, False)
# Initialize the aligner
aligner = Align.PairwiseAligner()
aligner.mode = 'global'  # For global alignment
aligner.match_score = alpha   # Positive value for matches
aligner.mismatch_score = -alpha  # Negative value for mismatches
aligner.open_gap_score = -beta  # Negative value for opening a gap
aligner.extend_gap_score = -beta
center_seq = sequences.pop(center)
center_seq, aligned_seqs= MSA(center_seq, sequences)

print("\nMutiple String Alignment")
print(f"{center}: {center_seq}")
for i in aligned_seqs.keys():
    print(f"{i}: {aligned_seqs[i]}")


