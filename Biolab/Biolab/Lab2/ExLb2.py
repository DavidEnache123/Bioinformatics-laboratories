import random

CODON_TABLE = {
    # U row
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    
    # C row
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    
    # A row
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    
    # G row
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}

def generate_random_rna(length: int = 30) -> str:
    """Generates a random sequence of G, C, U, A bases."""
    bases = ["G", "C", "U", "A"]
    # Ensure sequence length is divisible by 3
    length = length - (length % 3)
    return "".join(random.choice(bases) for _ in range(length))

def translate_rna(rna_sequence: str) -> list[str]:
    """Translates RNA codons into their corresponding amino acids."""
    amino_acids = []
    
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        amino_acid = CODON_TABLE.get(codon, "Unknown")
        
        # Stop translation immediately if a STOP codon is reached
        if amino_acid == "STOP":
            amino_acids.append("STOP")
            break
            
        amino_acids.append(amino_acid)
        
    return amino_acids


# 1. Generate and show the random sequence
random_rna = generate_random_rna(length=30)
print(f"Generated RNA Sequence:\n{random_rna}\n")

# 2. Translate and show the protein result
protein = translate_rna(random_rna)
print("Translated Protein Chain:")
if protein:
    print(" - ".join(protein))
else:
    print("(No amino acids produced)")