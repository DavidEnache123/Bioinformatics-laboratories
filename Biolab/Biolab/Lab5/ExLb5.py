import random

# DNA helper functions
def random_dna(length):
    return "".join(random.choices("ACGT", k=length))

def reverse_complement(seq):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[b] for b in reversed(seq))

def build_transposon(tir_len=6, body_len=20):
    """
    Constructs a transposable element:
    TIR_left + Internal_Body + TIR_right (reverse complement of TIR_left)
    """
    tir_left = random_dna(tir_len)
    tir_right = reverse_complement(tir_left)
    body = random_dna(body_len)
    return tir_left + body + tir_right, tir_len

def insert_transposon(host_seq, te_seq, tsd_len=4, insert_pos=None):
    """
    Simulates insertion creating Target Site Duplications (Direct Repeats):
    ...Host_Left... + [TSD + TE + TSD] + ...Host_Right...
    """
    if insert_pos is None:
        insert_pos = random.randint(0, len(host_seq) - tsd_len)
    
    target_site = host_seq[insert_pos : insert_pos + tsd_len]
    
    # Left flank up to target site
    left_part = host_seq[:insert_pos]
    # Right flank from after the target site
    right_part = host_seq[insert_pos + tsd_len:]
    
    # After insertion, target site is duplicated on both sides of the TE
    new_seq = left_part + target_site + te_seq + target_site + right_part
    te_start = insert_pos
    te_end = insert_pos + len(target_site) * 2 + len(te_seq)
    
    return new_seq, te_start, te_end, target_site

# Detection Algorithm

def detect_transposons(dna_seq, min_te_len=20, max_te_len=80, tsd_len=4, tir_len=6):
    """
    Scans the DNA sequence to locate:
    Direct Repeat (TSD) + Left TIR + ...Body... + Right TIR + Direct Repeat (TSD)
    """
    detected = []
    n = len(dna_seq)
    
    for i in range(n - (min_te_len + 2 * tsd_len)):
        tsd1 = dna_seq[i : i + tsd_len]
        
        # Scan possible transposon lengths
        for length in range(min_te_len, min(max_te_len, n - i - 2 * tsd_len) + 1):
            j = i + tsd_len + length
            tsd2 = dna_seq[j : j + tsd_len]
            
            # Check 1: Target Site Duplication (Direct Repeats match)
            if tsd1 == tsd2:
                te_candidate = dna_seq[i + tsd_len : j]
                tir_left = te_candidate[:tir_len]
                tir_right = te_candidate[-tir_len:]
                
                # Check 2: Terminal Inverted Repeats match (reverse complements)
                if tir_right == reverse_complement(tir_left):
                    detected.append({
                        "start": i,
                        "end": j + tsd_len,
                        "tsd": tsd1,
                        "tir": tir_left,
                        "te_length": len(te_candidate)
                    })
                    
    return detected


# Simulation & Experiment

if __name__ == "__main__":
    random.seed(42)  # For reproducible demonstration
    
    # 1. Generate base DNA host of ~250 bp
    host_dna = random_dna(250)
    print(f"Initial Host DNA Length: {len(host_dna)} bp\n")
    
    ground_truth = []
    current_seq = host_dna
    
    # Simulate 3-4 Transposable Elements
    num_elements = 3
    print("--- Simulating Insertions ---")
    
    # Element 1 & 2: Standard non-intersecting insertions
    for k in range(2):
        te, _ = build_transposon(tir_len=6, body_len=24)
        current_seq, start, end, tsd = insert_transposon(current_seq, te, tsd_len=4)
        ground_truth.append((f"TE_{k+1}", start, end, tsd))
        print(f"Inserted TE_{k+1} at index [{start}:{end}] with TSD '{tsd}'")

    # Element 3: Intersected / Nested insertion (inserted into an existing TE)
    # Target insertion inside TE_2
    nested_target = ground_truth[-1][1] + 12
    te_nested, _ = build_transposon(tir_len=6, body_len=18)
    current_seq, start, end, tsd = insert_transposon(current_seq, te_nested, tsd_len=4, insert_pos=nested_target)
    ground_truth.append(("TE_3 (Nested)", start, end, tsd))
    print(f"Inserted TE_3 (Nested inside TE_2) at index [{start}:{end}] with TSD '{tsd}'")
    
    print(f"\nFinal Artificial Sequence Length: {len(current_seq)} bp")
    print(f"Final Sequence (first 100 bp): {current_seq[:100]}...\n")
    
    # 2. Run Detection Software
    print("--- Running Detection Application ---")
    detections = detect_transposons(current_seq, min_te_len=20, max_te_len=80, tsd_len=4, tir_len=6)
    
    print(f"Detected {len(detections)} candidate element(s):")
    for idx, d in enumerate(detections, 1):
        print(f"  {idx}. Start: {d['start']}, End: {d['end']} (Length: {d['end'] - d['start']} bp) | TSD: {d['tsd']} | Left TIR: {d['tir']}")