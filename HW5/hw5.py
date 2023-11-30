def read_fasta(file_path):
    """
    Read sequences from a FASTA file.
    """
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

def cost_function(char1, char2, alpha, beta):
    """
    Calculate the cost of aligning two characters, based on the given cost function parameters.
    """
    if char1 == char2:
        return 0
    elif char1 == '-' or char2 == '-':
        return beta
    else:
        return alpha

def pairwise_distance(seq1, seq2, alpha, beta):
    """
    Calculate the pairwise distance between two sequences, based on the given cost function parameters.
    """
    # Pad the shorter sequence with gaps if the lengths are unequal
    l1, l2 = len(seq1), len(seq2)
    if l1 < l2:
        seq1 += '-' * (l2 - l1)
    elif l2 < l1:
        seq2 += '-' * (l1 - l2)

    distance = 0
    for a, b in zip(seq1, seq2):
        distance += cost_function(a, b, alpha, beta)
    return distance

def find_center_sequence(sequences, alpha, beta):
    """
    Find the center sequence from a list of sequences, based on the total pairwise distances.
    """
    min_distance = float('inf')
    center_sequence = None
    for seq in sequences:
        total_distance = sum(pairwise_distance(seq, other, alpha, beta)
                             for other in sequences if other != seq)
        if total_distance < min_distance:
            min_distance = total_distance
            center_sequence = seq
    return center_sequence

# Example sequences
sequences = ["acctt", "tcggc", "tactt", "atcgt", "acata"]
alpha = 1  # Cost for mismatch
beta = 2   # Cost for gap

# Find the center sequence
center_sequence = find_center_sequence(sequences, alpha, beta)
print(center_sequence)

