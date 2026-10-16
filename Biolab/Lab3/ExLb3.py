from collections import Counter
from pathlib import Path
import sys

# Standard Genetic Code Dictionary (RNA codon -> Amino Acid)
CODON_TABLE = {
    "UUU": "Phenylalanine (Phe)",
    "UUC": "Phenylalanine (Phe)",
    "UUA": "Leucine (Leu)",
    "UUG": "Leucine (Leu)",
    "CUU": "Leucine (Leu)",
    "CUC": "Leucine (Leu)",
    "CUA": "Leucine (Leu)",
    "CUG": "Leucine (Leu)",
    "AUU": "Isoleucine (Ile)",
    "AUC": "Isoleucine (Ile)",
    "AUA": "Isoleucine (Ile)",
    "AUG": "Methionine (Met)",
    "GUU": "Valine (Val)",
    "GUC": "Valine (Val)",
    "GUA": "Valine (Val)",
    "GUG": "Valine (Val)",
    "UCU": "Serine (Ser)",
    "UCC": "Serine (Ser)",
    "UCA": "Serine (Ser)",
    "UCG": "Serine (Ser)",
    "CCU": "Proline (Pro)",
    "CCC": "Proline (Pro)",
    "CCA": "Proline (Pro)",
    "CCG": "Proline (Pro)",
    "ACU": "Threonine (Thr)",
    "ACC": "Threonine (Thr)",
    "ACA": "Threonine (Thr)",
    "ACG": "Threonine (Thr)",
    "GCU": "Alanine (Ala)",
    "GCC": "Alanine (Ala)",
    "GCA": "Alanine (Ala)",
    "GCG": "Alanine (Ala)",
    "UAU": "Tyrosine (Tyr)",
    "UAC": "Tyrosine (Tyr)",
    "UAA": "Stop",
    "UAG": "Stop",
    "UGA": "Stop",
    "CAU": "Histidine (His)",
    "CAC": "Histidine (His)",
    "CAA": "Glutamine (Gln)",
    "CAG": "Glutamine (Gln)",
    "AAU": "Asparagine (Asn)",
    "AAC": "Asparagine (Asn)",
    "AAA": "Lysine (Lys)",
    "AAG": "Lysine (Lys)",
    "GAU": "Aspartic acid (Asp)",
    "GAC": "Aspartic acid (Asp)",
    "GAA": "Glutamic acid (Glu)",
    "GAG": "Glutamic acid (Glu)",
    "UGU": "Cysteine (Cys)",
    "UGC": "Cysteine (Cys)",
    "UGG": "Tryptophan (Trp)",
    "CGU": "Arginine (Arg)",
    "CGC": "Arginine (Arg)",
    "CGA": "Arginine (Arg)",
    "CGG": "Arginine (Arg)",
    "AGU": "Serine (Ser)",
    "AGC": "Serine (Ser)",
    "AGA": "Arginine (Arg)",
    "AGG": "Arginine (Arg)",
    "GGU": "Glycine (Gly)",
    "GGC": "Glycine (Gly)",
    "GGA": "Glycine (Gly)",
    "GGG": "Glycine (Gly)",
}

FOOD_SOURCES = {
    "Leucine (Leu)": "Eggs, beef, chicken breast, tuna, soy, parmesan cheese.",
    "Isoleucine (Ile)": "Lentils, spirulina, turkey, eggs, cottage cheese, pumpkin seeds.",
    "Valine (Val)": "Mushrooms, peanuts, soy protein, beef, cheese, oats.",
    "Methionine (Met)": "Brazil nuts, eggs, fish, sesame seeds, poultry, dairy.",
    "Phenylalanine (Phe)": "Beef, poultry, pork, dairy, soybeans, pumpkin seeds.",
    "Threonine (Thr)": "Cottage cheese, lentils, sesame seeds, lean meats, poultry.",
    "Tryptophan (Trp)": "Turkey, whole milk, oats, chocolate, tofu, chia seeds.",
    "Lysine (Lys)": "Red meat, pork, parmesan cheese, cod, eggs, legumes, fenugreek.",
    "Histidine (His)": "Meat, whole wheat, rye, fish, beans, eggs.",
    "Alanine (Ala)": "Meat, fish, eggs, dairy, whole grains, avocados, beans.",
    "Arginine (Arg)": "Turkey, pumpkin seeds, peanuts, soybeans, dairy, chickpeas.",
    "Asparagine (Asn)": "Dairy, asparagus, potatoes, beef, legumes, whole grains.",
    "Aspartic acid (Asp)": "Sprouted seeds, beef, asparagus, poultry, avocado, molasses.",
    "Cysteine (Cys)": "Garlic, onions, egg yolks, oats, red peppers, poultry.",
    "Glutamic acid (Glu)": "Tomatoes, aged cheese, soy sauce, mushrooms, cured meats.",
    "Glutamine (Gln)": "Bone broth, cabbage, spinach, seafood, beef, ricotta.",
    "Glycine (Gly)": "Bone broth, gelatin, pork skin, chicken skin, legumes.",
    "Proline (Pro)": "Gelatin, bone broth, cottage cheese, cabbage, asparagus.",
    "Serine (Ser)": "Soybeans, eggs, peanuts, lentils, shellfish, dairy.",
    "Tyrosine (Tyr)": "Soy products, turkey, cheese, almonds, bananas, avocados.",
    "Stop": "N/A (Terminates translation; not an amino acid)",
}


def read_fasta(filepath: Path) -> str:
    sequence_lines = []
    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_lines.append(line.upper())
    return "".join(sequence_lines)


def get_codon_frequencies(sequence: str, step: int = 3) -> tuple[Counter, int]:
    codon_counts = Counter()
    total_codons = 0

    for i in range(0, len(sequence) - 2, step):
        codon = sequence[i : i + 3]
        if len(codon) == 3 and all(base in "ATCG" for base in codon):
            codon_counts[codon] += 1
            total_codons += 1

    return codon_counts, total_codons


def translate_dna_to_rna(dna_codon: str) -> str:
    return dna_codon.replace("T", "U")


def main():
    # Automatically locate sequence.fasta in the exact same directory as this script
    script_dir = Path(__file__).resolve().parent
    fasta_file = script_dir / "sequence.fasta"

    if not fasta_file.exists():
        print(f"Error: Could not find '{fasta_file.name}' in {script_dir}")
        sys.exit(1)

    print(f"Reading: {fasta_file}")
    dna_sequence = read_fasta(fasta_file)

    if not dna_sequence:
        print("Error: The FASTA file contains no sequence data.")
        sys.exit(1)

    print(f"Loaded genome with {len(dna_sequence):,} base pairs.")

    codon_counts, total_valid_codons = get_codon_frequencies(
        dna_sequence, step=3
    )

    if total_valid_codons == 0:
        print("No valid 3-letter codons found in the sequence.")
        sys.exit(1)

    top_10 = codon_counts.most_common(10)

    print(f"\nTotal codons processed: {total_valid_codons:,}")
    print(f"{'Rank':<5} {'DNA':<6} {'RNA':<6} {'Count':<10} {'Freq (%)':<10}")
    print("-" * 40)

    for rank, (dna_codon, count) in enumerate(top_10, 1):
        rna_codon = translate_dna_to_rna(dna_codon)
        freq = (count / total_valid_codons) * 100
        print(
            f"{rank:<5} {dna_codon:<6} {rna_codon:<6} {count:<10} {freq:<10.2f}"
        )

    print("\n" + "=" * 80)
    print("TOP 10 AMINO ACIDS & DIETARY NUTRIENT SOURCES")
    print("=" * 80)

    seen_aa = set()
    for rank, (dna_codon, _) in enumerate(top_10, 1):
        rna_codon = translate_dna_to_rna(dna_codon)
        amino_acid = CODON_TABLE.get(rna_codon, "Unknown")
        sources = FOOD_SOURCES.get(
            amino_acid, "Dietary data currently unavailable."
        )

        print(
            f"#{rank} Codon: {dna_codon} -> RNA: {rna_codon} => {amino_acid}"
        )
        if amino_acid not in seen_aa:
            print(f"   Dietary Sources: {sources}")
            seen_aa.add(amino_acid)
        else:
            print("   Dietary Sources: (Same as previously listed above)")
        print("-" * 80)


if __name__ == "__main__":
    main()