from pathlib import Path
import random
import sys


def read_fasta(filepath: Path) -> str:
    """Reads a FASTA file and returns the concatenated nucleotide sequence in uppercase."""
    sequence_lines = []
    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_lines.append(line.upper())
    return "".join(sequence_lines)


def sample_genome_uniformly(
    sequence: str,
    num_samples: int = 20_000,
    min_len: int = 100,
    max_len: int = 150,
) -> list[tuple[int, int, str]]:
    """Samples fragments uniformly across the genome.

    Returns a list of tuples: (start_index_1_based, length, sequence_slice)
    """
    seq_length = len(sequence)

    if seq_length < max_len:
        raise ValueError(
            f"Genome length ({seq_length}) is smaller than maximum sample length ({max_len})."
        )

    samples = []

    for _ in range(num_samples):
        # 1. Random sample length uniformly between 100 and 150
        sample_len = random.randint(min_len, max_len)

        # 2. Random starting index chosen uniformly across all valid start positions
        max_start = seq_length - sample_len
        start_idx = random.randint(0, max_start)

        # 3. Slice the fragment
        fragment = sequence[start_idx : start_idx + sample_len]

        # Storing with 1-based start index for standard genomic convention
        samples.append((start_idx + 1, sample_len, fragment))

    return samples


def save_samples_to_fasta(
    samples: list[tuple[int, int, str]], output_path: Path
):
    """Saves generated genomic samples to a new FASTA file."""
    with open(output_path, "w") as out:
        for i, (start_pos, length, frag) in enumerate(samples, 1):
            # FASTA header: >sample_ID start=START length=LEN
            out.write(f">sample_{i} start_pos={start_pos} length={length}\n")
            out.write(f"{frag}\n")


def main():
    # Resolve the path to sequence.fasta in the current script's folder (Lab4)
    script_dir = Path(__file__).resolve().parent
    fasta_file = script_dir / "sequence.fasta"
    output_file = script_dir / "samples_20k.fasta"

    if not fasta_file.exists():
        print(f"Error: Could not find '{fasta_file.name}' in {script_dir}")
        sys.exit(1)

    print(f"Reading genome from: {fasta_file.name}...")
    dna_sequence = read_fasta(fasta_file)
    genome_len = len(dna_sequence)
    print(f"Loaded genome with {genome_len:,} base pairs.")

    num_samples = 20_000
    min_len = 100
    max_len = 150

    print(
        f"Generating {num_samples:,} uniform samples (length: {min_len}-{max_len} bp)..."
    )
    samples = sample_genome_uniformly(
        dna_sequence,
        num_samples=num_samples,
        min_len=min_len,
        max_len=max_len,
    )

    print(f"Saving samples to: {output_file.name}...")
    save_samples_to_fasta(samples, output_file)

    print(f"\nSampling Complete!")
    print(f"- Total samples collected: {len(samples):,}")
    print(f"- Output saved to: {output_file.resolve()}")

    # Display the first 3 samples as a preview
    print("\n" + "=" * 60)
    print("PREVIEW OF FIRST 3 SAMPLES:")
    print("=" * 60)
    for i in range(min(3, len(samples))):
        start_pos, length, frag = samples[i]
        print(f"Sample #{i+1} | Start: {start_pos:,} | Length: {length} bp")
        print(f"Sequence: {frag[:60]}... ({length} bp total)\n")


if __name__ == "__main__":
    main()