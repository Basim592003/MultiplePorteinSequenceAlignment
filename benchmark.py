import numpy as np
from Bio import SeqIO

def calculate_percent_identity(seq1, seq2):
    """Calculate the percent identity between two aligned sequences, ignoring gaps."""
    matches, total = 0, 0

    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':  
            total += 1
            if a == b:
                matches += 1

    return (matches / total * 100) if total > 0 else 0.0

def create_identity_matrix(records):
    """Create a percent identity matrix from sequence records."""
    num_sequences = len(records)
    matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        seq1 = str(records[i].seq)
        for j in range(i, num_sequences):
            seq2 = str(records[j].seq)
            percent_identity = calculate_percent_identity(seq1, seq2)
            matrix[i][j] = matrix[j][i] = round(percent_identity, 2) 

    return matrix

def evaluate_alignment(algorithm, alignment_file):
    """Evaluate the alignment results, returning various metrics."""
    records = list(SeqIO.parse(alignment_file, "fasta"))

    total_sequences = len(records)
    total_length = max(len(record.seq) for record in records)
    gap_count = sum(str(record.seq).count('-') for record in records)

    identity_matrix = create_identity_matrix(records)
    sequence_names = [record.id for record in records]  # Get sequence names

    return {
        'algorithm': algorithm,
        'total_sequences': total_sequences,
        'total_length': total_length,
        'gap_count': gap_count,
        'identity_matrix': identity_matrix.tolist(),
        'sequence_names': sequence_names  # Add sequence names to the return dict
    }
