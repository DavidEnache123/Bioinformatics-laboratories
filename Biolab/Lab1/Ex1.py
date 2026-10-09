from collections import Counter


def analyze_sequence_kmers(
    sequence: str, k: int
) -> list[tuple[str, int, float]]:
    """Calculates counts and percentage frequencies for any combination length (k).

    Returns a list of (pattern, count, percentage) sorted in descending order.
    """
    clean_seq = "".join(sequence.upper().split())
    n = len(clean_seq)
    total_windows = n - k + 1

    if total_windows <= 0:
        return []

    # Extract all overlapping substrings of length k
    kmers = [clean_seq[i : i + k] for i in range(total_windows)]
    counts = Counter(kmers)

    # Calculate frequencies and sort: highest count first, then alphabetically
    results = [
        (kmer, count, (count / total_windows) * 100)
        for kmer, count in counts.items()
    ]
    results.sort(key=lambda item: (-item[1], item[0]))

    return results


def display_chart(
    data: list[tuple[str, int, float]],
    title: str,
    unit_label: str = "Pattern",
):
    """Renders a formatted table and an ASCII horizontal bar chart."""
    print(f"\n{'=' * 56}")
    print(f"{title}")
    print(f"{'=' * 56}")

    if not data:
        print("Sequence too short for this length.")
        return

    print(
        f"{unit_label:<8} | {'Count':<6} | {'Freq (%)':<9} | {'Distribution Chart'}"
    )
    print("-" * 56)

    max_bar_width = 20
    max_freq = data[0][2] if data else 1

    for pattern, count, freq in data:
        bar_len = int((freq / max_freq) * max_bar_width) if max_freq > 0 else 0
        bar = "█" * max(bar_len, 1 if count > 0 else 0)
        print(f"{pattern:<8} | {count:<6} | {freq:>7.2f}% | {bar}")


def run_full_dna_analysis(sequence: str):
    clean_seq = "".join(sequence.upper().split())
    print(f"Sequence: {clean_seq}")
    print(f"Total Length: {len(clean_seq)} characters")

    # 1. Single letter frequencies (k = 1)
    single_letters = analyze_sequence_kmers(clean_seq, k=1)
    display_chart(
        single_letters,
        "1-Letter Frequencies (Descending)",
        unit_label="Letter",
    )

    # 2. Dinucleotides / 2-letter combinations (k = 2)
    two_letters = analyze_sequence_kmers(clean_seq, k=2)
    display_chart(
        two_letters, "2-Letter Combinations (Descending)", unit_label="Dimer"
    )

    # 3. Trinucleotides / 3-letter combinations (k = 3)
    three_letters = analyze_sequence_kmers(clean_seq, k=3)
    display_chart(
        three_letters, "3-Letter Combinations (Descending)", unit_label="Trimer"
    )


if __name__ == "__main__":
    # Test sequence containing standard bases and the 'EE' test cases
    dna_input = "AGCTTTTCATTIITGACTGCA"
    run_full_dna_analysis(dna_input)