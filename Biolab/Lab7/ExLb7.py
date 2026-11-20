import os
import random


def load_fasta_sequence(filepath):
    """Reads the full sequence from a single-record FASTA file."""
    header = ""
    seq_chunks = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
            else:
                seq_chunks.append(line.upper())

    return header, "".join(seq_chunks)


def get_random_subsequences(genome, min_len=10, max_len=20):
    """Extracts two random subsequences of random lengths from the genome."""
    genome_len = len(genome)

    # Pick random lengths for S1 and S2
    len1 = random.randint(min_len, max_len)
    len2 = random.randint(min_len, max_len)

    # Pick random starting indices
    start1 = random.randint(0, genome_len - len1)
    start2 = random.randint(0, genome_len - len2)

    s1 = genome[start1 : start1 + len1]
    s2 = genome[start2 : start2 + len2]

    return (start1, s1), (start2, s2)


def global_alignment(s1, s2, match=2, mismatch=-1, gap=-2):
    """Performs Needleman-Wunsch global sequence alignment."""
    n = len(s2)  # Rows
    m = len(s1)  # Columns

    # 1. Initialize scoring matrix
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # 2. Fill first column and first row with gap penalties
    for i in range(n + 1):
        score_matrix[i][0] = i * gap
    for j in range(m + 1):
        score_matrix[0][j] = j * gap

    # 3. Fill the rest of the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score_matrix[i - 1][j - 1] + (match if s2[i - 1] == s1[j - 1] else mismatch)
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(diag, up, left)

    # 4. Traceback from bottom-right to top-left
    aligned_s1 = ""
    aligned_s2 = ""
    i, j = n, m

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]

        if i > 0 and j > 0:
            step_score = match if s2[i - 1] == s1[j - 1] else mismatch
            if current_score == score_matrix[i - 1][j - 1] + step_score:
                aligned_s1 = s1[j - 1] + aligned_s1
                aligned_s2 = s2[i - 1] + aligned_s2
                i -= 1
                j -= 1
                continue

        if i > 0 and current_score == score_matrix[i - 1][j] + gap:
            aligned_s1 = "-" + aligned_s1
            aligned_s2 = s2[i - 1] + aligned_s2
            i -= 1
        else:
            aligned_s1 = s1[j - 1] + aligned_s1
            aligned_s2 = "-" + aligned_s2
            j -= 1

    return score_matrix, aligned_s1, aligned_s2


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_path = os.path.join(script_dir, "sequence.fasta")

    # Load the genome and extract two random subsequences
    header, genome = load_fasta_sequence(fasta_path)
    (pos1, S1), (pos2, S2) = get_random_subsequences(genome, min_len=8, max_len=16)

    print(f"Source: {header}")
    print(f"S1 (index {pos1:,} to {pos1 + len(S1):,}): {S1}")
    print(f"S2 (index {pos2:,} to {pos2 + len(S2):,}): {S2}")

    # Run Needleman-Wunsch with the whiteboard scoring scheme
    matrix, align1, align2 = global_alignment(S1, S2, match=2, mismatch=-1, gap=-2)

    # Print the Scoring Matrix
    print("\nScoring Matrix:")
    col_header = [" ", "-"] + list(S1)
    print("\t".join(col_header))
    for label, row in zip(["-"] + list(S2), matrix):
        print(f"{label}\t" + "\t".join(map(str, row)))

    # Print Alignment Result
    match_line = "".join("|" if a == b else " " for a, b in zip(align1, align2))
    print("\nOptimal Global Alignment:")
    print(f"S1: {align1}")
    print(f"    {match_line}")
    print(f"S2: {align2}")
    print(f"\nOptimal Alignment Score: {matrix[-1][-1]}")