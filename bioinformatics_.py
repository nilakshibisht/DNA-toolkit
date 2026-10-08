

#!/usr/bin/env python3
"""
DNA Toolkit - a simple bioinformatics starter project.
 
Features:
  - Validate a DNA sequence
  - Count nucleotides
  - GC content
  - Complement and reverse complement
  - Transcription (DNA -> RNA)
  - Translation (DNA -> protein)
 
Run:  python3 dna_toolkit.py
"""
 
# Standard codon table (RNA codons). '*' means STOP.
CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
 
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
 
 
def clean_sequence(seq: str) -> str:
    """Remove whitespace/newlines and convert to uppercase."""
    return "".join(seq.split()).upper()
 
 
def is_valid_dna(seq: str) -> bool:
    """Return True if the sequence contains only A, T, G, C."""
    return len(seq) > 0 and set(seq) <= set("ATGC")
 
 
def count_nucleotides(seq: str) -> dict:
    """Count each nucleotide."""
    return {base: seq.count(base) for base in "ATGC"}
 
 
def gc_content(seq: str) -> float:
    """Return GC content as a percentage."""
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 2)
 
 
def complement(seq: str) -> str:
    """Return the complement strand (same direction)."""
    return "".join(COMPLEMENT[base] for base in seq)
 
 
def reverse_complement(seq: str) -> str:
    """Return the reverse complement (the opposite strand, read 5' -> 3')."""
    return complement(seq)[::-1]
 
 
def transcribe(seq: str) -> str:
    """Transcribe DNA (coding strand) into mRNA by replacing T with U."""
    return seq.replace("T", "U")
 
 
def translate(seq: str, stop_at_stop_codon: bool = True) -> str:
    """Translate a DNA sequence into a protein sequence."""
    rna = transcribe(seq)
    protein = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i + 3]
        amino_acid = CODON_TABLE[codon]
        if amino_acid == "*" and stop_at_stop_codon:
            break
        protein.append(amino_acid)
    return "".join(protein)
 
 
def analyze(seq: str) -> None:
    """Run all analyses and print a report."""
    seq = clean_sequence(seq)
 
    if not is_valid_dna(seq):
        print("Error: sequence must contain only A, T, G, C.")
        return
 
    print("\n" + "=" * 50)
    print("DNA ANALYSIS REPORT")
    print("=" * 50)
    print(f"Sequence        : {seq}")
    print(f"Length          : {len(seq)} bp")
    print(f"Nucleotides     : {count_nucleotides(seq)}")
    print(f"GC content      : {gc_content(seq)} %")
    print(f"Complement      : {complement(seq)}")
    print(f"Reverse compl.  : {reverse_complement(seq)}")
    print(f"mRNA            : {transcribe(seq)}")
    print(f"Protein         : {translate(seq)}")
    if len(seq) % 3 != 0:
        print("Note: length is not a multiple of 3; trailing bases were ignored in translation.")
    print("=" * 50 + "\n")
 
 
def main() -> None:
    print("DNA Toolkit (type 'q' to quit)")
    print("Example: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
 
    while True:
        user_input = input("\nEnter a DNA sequence: ")
        if user_input.strip().lower() == "q":
            print("Goodbye!")
            break
        analyze(user_input)
 
 
if __name__ == "__main__":
    main()
 

