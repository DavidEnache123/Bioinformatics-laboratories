import os
import sys
from collections import defaultdict

# 1. Load Local FASTA Sequence

def load_fasta_sequence(filepath):
    """Parses a FASTA file and returns the concatenated uppercase sequence."""
    if not os.path.exists(filepath):
        print(f"Error: Could not find '{filepath}'. Check the file path.")
        sys.exit(1)

    seq_chunks = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_chunks.append(line.upper())
    return "".join(seq_chunks)


# 2. Sequence Utilities

COMPLEMENT_MAP = str.maketrans("ACGTN", "TGCAN")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT_MAP)[::-1]


# 3. Two-Window Sliding Algorithm (S1 & S2)

def find_transposons(
    seq,
    window_len=12,         # Length of S1 and S2 (unknown TIR size)
    min_distance=100,      # Minimum span between S1 and S2
    max_distance=15000,    # Max internal distance to keep search focused
    max_matches=50         # Stop early once enough candidates are found
):
    """
    S1 slides across the genome step-by-step.
    S2 matches the reverse complement of S1 downstream, separated by >= 100 bp.
    """
    n = len(seq)
    print(f"Loaded sequence length: {n:,} bp")
    print(f"Indexing k-mers (length={window_len})...")

    # Index position of all k-mers to keep the S2 search near O(1)
    kmer_positions = defaultdict(list)
    for pos in range(n - window_len + 1):
        kmer = seq[pos : pos + window_len]
        if "N" not in kmer:
            kmer_positions[kmer].append(pos)

    print("Sliding S1 and scanning for inverted complementary S2 downstream...")
    results = []

    # Slide S1 across the sequence with step = 1
    for s1_start in range(n - window_len - min_distance):
        s1_seq = seq[s1_start : s1_start + window_len]
        if "N" in s1_seq:
            continue

        # Inverted and complementary sequence that S2 must match
        target_s2 = reverse_complement(s1_seq)

        if target_s2 in kmer_positions:
            for s2_start in kmer_positions[target_s2]:
                dist = s2_start - (s1_start + window_len)
                if min_distance <= dist <= max_distance:
                    s2_end = s2_start + window_len
                    total_span = s2_end - s1_start

                    results.append({
                        "s1_start": s1_start,
                        "s1_end": s1_start + window_len,
                        "s1_seq": s1_seq,
                        "s2_start": s2_start,
                        "s2_end": s2_end,
                        "s2_seq": target_s2,
                        "distance_between": dist,
                        "total_span": total_span
                    })

                    if len(results) >= max_matches:
                        print(f"Reached cap of {max_matches} candidates. Stopping search.")
                        return results

    return results


# 4. Execution

if __name__ == "__main__":
    # Resolve the path relative to this script so it runs from anywhere
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_path = os.path.join(script_dir, "sequence.fasta")

    # Load sequence
    genome_seq = load_fasta_sequence(fasta_path)

    # Detect transposons (S1/S2 window = 12 bp, minimum distance = 100 bp)
    candidates = find_transposons(
        genome_seq,
        window_len=12,
        min_distance=100,
        max_distance=15000,
        max_matches=25
    )

    print(f"\nFound {len(candidates)} candidate inverted repeat regions:\n")
    print(f"{'#':<3} | {'S1 Range':<16} | {'S1 Seq':<14} | {'S2 Range':<16} | {'S2 Seq (RevComp)':<16} | {'Total Span':<10}")
    print("-" * 88)
    for idx, c in enumerate(candidates, 1):
        s1_range = f"[{c['s1_start']}:{c['s1_end']}]"
        s2_range = f"[{c['s2_start']}:{c['s2_end']}]"
        print(f"{idx:<3} | {s1_range:<16} | {c['s1_seq']:<14} | {s2_range:<16} | {c['s2_seq']:<16} | {c['total_span']:<10} bp")