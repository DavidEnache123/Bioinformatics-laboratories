from collections import defaultdict
from pathlib import Path
import sys
import time


def load_samples(filepath: Path) -> list[str]:
    """Reads the sampled FASTA file and strips out all metadata,

    retaining only the raw sample strings without access to the original genome.
    """
    reads = []
    current_seq = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq:
                    reads.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line.upper())
        if current_seq:
            reads.append("".join(current_seq))
    return reads


def assemble_greedy_overlap(
    reads: list[str], min_overlap: int = 30, kmer_size: int = 25
) -> list[str]:
    """Assembles overlapping read fragments into contigs using a greedy overlap layout.

    An indexed k-mer map is used to find candidate suffix-prefix matches efficiently.
    """
    # Remove duplicates and exact sub-reads to eliminate redundancy
    unique_reads = list(set(reads))
    unique_reads.sort(key=len, reverse=True)

    contigs = []
    for r in unique_reads:
        if not any(r in c for c in contigs):
            contigs.append(r)

    print(f"Starting assembly from {len(contigs):,} unique contigs...")

    # Build k-mer prefix index: maps a k-mer prefix to the list of contigs starting with it
    prefix_index = defaultdict(list)
    for i, c in enumerate(contigs):
        prefix_index[c[:kmer_size]].append(i)

    merged = True
    iteration = 0

    while merged:
        merged = False
        iteration += 1

        best_overlap = 0
        best_pair = (-1, -1)

        # Check contigs for suffix-to-prefix matches
        for i, c1 in enumerate(contigs):
            if not c1:
                continue

            # Look at candidate prefixes that match the end segments of c1
            suffix_kmer = c1[-kmer_size:]
            candidates = prefix_index.get(suffix_kmer, [])

            for j in candidates:
                if i == j:
                    continue
                c2 = contigs[j]
                if not c2:
                    continue

                # Find maximum suffix-prefix overlap between c1 and c2
                max_possible = min(len(c1), len(c2))
                for ov in range(
                    max_possible, max(min_overlap, best_overlap), -1
                ):
                    if c1.endswith(c2[:ov]):
                        best_overlap = ov
                        best_pair = (i, j)
                        break

            if best_overlap >= 80:
                # Greedy choice: accept a sufficiently long overlap immediately to speed up assembly
                break

        # Merge the best pair found
        if best_pair != (-1, -1):
            i, j = best_pair
            c1, c2 = contigs[i], contigs[j]
            merged_contig = c1 + c2[best_overlap:]

            # Overwrite i with merged contig, invalidate j
            contigs[i] = merged_contig
            contigs[j] = ""

            # Rebuild prefix index for active contigs
            prefix_index = defaultdict(list)
            for idx, c in enumerate(contigs):
                if c and len(c) >= kmer_size:
                    prefix_index[c[:kmer_size]].append(idx)

            merged = True
            if iteration % 200 == 0:
                active_count = sum(1 for c in contigs if c)
                print(
                    f"Iteration {iteration}: {active_count:,} contigs remaining..."
                )

    # Collect final active contigs, sorted by length descending
    assembled_contigs = [c for c in contigs if c]
    assembled_contigs.sort(key=len, reverse=True)
    return assembled_contigs


def save_assembled_contigs(contigs: list[str], output_path: Path):
    """Saves the reconstructed sequence(s) to a standard FASTA file."""
    with open(output_path, "w") as out:
        for idx, contig in enumerate(contigs, 1):
            out.write(f">assembled_contig_{idx} length={len(contig)} bp\n")
            for pos in range(0, len(contig), 80):
                out.write(contig[pos : pos + 80] + "\n")


def main():
    # Both ExLb4Pt2.py and samples_20k.fasta are directly in Lab4
    lab_dir = Path(__file__).resolve().parent
    samples_file = lab_dir / "samples_20k.fasta"
    output_file = lab_dir / "reconstructed_genome.fasta"

    if not samples_file.exists():
        print(f"Error: Could not find '{samples_file.name}' in {lab_dir}")
        sys.exit(1)

    print(f"Loading samples from: {samples_file.name}...")
    start_time = time.time()
    reads = load_samples(samples_file)
    print(f"Loaded {len(reads):,} sequence fragments.")

    print("\nReconstructing genome via de novo greedy overlap assembly...")
    contigs = assemble_greedy_overlap(reads, min_overlap=30, kmer_size=20)

    elapsed = time.time() - start_time
    print(f"\nAssembly finished in {elapsed:.2f} seconds.")
    print(f"Total reconstructed contigs: {len(contigs):,}")

    if contigs:
        primary_contig = contigs[0]
        print(f"Longest assembled contig: {len(primary_contig):,} bp")
        save_assembled_contigs(contigs, output_file)
        print(f"Saved reconstructed genome to: {output_file.name}")
    else:
        print("Could not assemble contigs with the given overlap criteria.")


if __name__ == "__main__":
    main()