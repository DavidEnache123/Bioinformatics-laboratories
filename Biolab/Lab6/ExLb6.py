import os
import re
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# 1. Define the 5 restriction enzymes: (Recognition Sequence, Cut Offset on 5'->3' strand)
ENZYMES = {
    "EcoRI":   ("GAATTC", 1),  # G^AATTC
    "BamHI":   ("GGATCC", 1),  # G^GATCC
    "HindIII": ("AAGCTT", 1),  # A^AGCTT
    "TaqI":    ("TCGA",   1),  # T^CGA
    "HaeIII":  ("GGCC",   2),  # GG^CC
}

def load_dna_sequence(fasta_path):
    """Reads a FASTA file from NCBI and returns a clean uppercase DNA string."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Could not find '{fasta_path}'. Make sure sequence.fasta is in the Lab6 folder.")
    
    with open(fasta_path, "r") as f:
        header = f.readline().strip()
        print(f"Loaded Sequence: {header}")
        lines = [line.strip() for line in f if not line.startswith(">")]
    return "".join(lines).upper()

def digest_dna(dna_seq, enzymes):
    """Finds cleavage sites and calculates fragment lengths for each enzyme."""
    results = {}
    seq_len = len(dna_seq)

    for name, (site, offset) in enzymes.items():
        cut_positions = [
            match.start() + offset 
            for match in re.finditer(f"(?={site})", dna_seq)
        ]
        
        boundaries = [0] + cut_positions + [seq_len]
        fragments = [
            boundaries[i + 1] - boundaries[i] 
            for i in range(len(boundaries) - 1)
        ]
        
        results[name] = {
            "cuts": len(cut_positions),
            "positions": cut_positions,
            "fragments": sorted(fragments, reverse=True)
        }
    return results

def format_preview(items, limit=10):
    """Formats long lists cleanly for terminal output without flooding the console."""
    if len(items) <= limit:
        return str(items)
    preview = ", ".join(map(str, items[:limit]))
    return f"[{preview}, ... ({len(items) - limit} more)]"

def simulate_gel(results):
    """Simulates an agarose gel electrophoresis plot optimized for full bacterial genomes."""
    fig, ax = plt.subplots(figsize=(10, 6.5))
    
    # Dark gel background styling
    ax.set_facecolor("#121212")
    fig.patch.set_facecolor("#1e1e1e")

    # Extended DNA Ladder suitable for full genomic digests (in bp)
    ladder = [50000, 25000, 10000, 5000, 2500, 1000, 500, 250, 100, 50]
    lanes = ["Ladder"] + list(results.keys())

    # Plot Ladder bands (Lane 0)
    for size in ladder:
        ax.hlines(y=size, xmin=-0.3, xmax=0.3, color="#00ffcc", linewidth=2.5, alpha=0.9)
        ax.text(-0.42, size, f"{size:,} bp", color="white", va="center", ha="right", fontsize=8)

    # Plot Enzyme Digestion bands using LineCollection for fast rendering of thousands of bands
    all_segments = []
    max_frag = max(ladder)

    for idx, (enzyme, data) in enumerate(results.items(), start=1):
        for frag in data["fragments"]:
            if frag > 0:
                all_segments.append([(idx - 0.3, frag), (idx + 0.3, frag)])
                if frag > max_frag:
                    max_frag = frag

    # Semi-transparent green bands simulate fluorescence intensity (dense regions form smears)
    band_collection = LineCollection(all_segments, colors="#39ff14", linewidths=1.5, alpha=0.15)
    ax.add_collection(band_collection)

    # Configure logarithmic Y-axis
    ax.set_yscale("log")
    ax.set_ylim(30, max_frag * 1.4)
    ax.set_xlim(-1.1, len(lanes) - 0.3)
    
    # Formatting axes and labels
    ax.set_xticks(range(len(lanes)))
    ax.set_xticklabels(lanes, color="white", fontsize=11, fontweight="bold")
    ax.xaxis.tick_top()  # Wells are at the top of the gel
    ax.set_ylabel("Fragment Size (bp) - Log Scale", color="white", fontsize=11)
    ax.set_title("Restriction Enzyme Digestion - E. coli K-12 Gel Simulation", color="white", pad=20, fontsize=13)
    
    ax.tick_params(axis="y", which="both", left=False, labelleft=False)
    for spine in ax.spines.values():
        spine.set_edgecolor("#555555")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    fasta_file = os.path.join(os.path.dirname(__file__), "sequence.fasta")
    dna = load_dna_sequence(fasta_file)
    
    print(f"Total Genome Length: {len(dna):,} bp\n" + "-" * 60)
    
    digestion_results = digest_dna(dna, ENZYMES)
    
    for enzyme, info in digestion_results.items():
        print(f"Enzyme: {enzyme} ({ENZYMES[enzyme][0]})")
        print(f"  1. Number of cleavages     : {info['cuts']:,}")
        print(f"  2. Cleavage positions      : {format_preview(info['positions'])}")
        print(f"  3. Largest fragments (bp)  : {format_preview(info['fragments'])}\n")
        
    simulate_gel(digestion_results)