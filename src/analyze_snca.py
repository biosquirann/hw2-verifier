#!/usr/bin/env python3
"""Analyze SNCA protein: length and amino acid composition."""

from pathlib import Path

FASTA_PATH = Path(__file__).resolve().parent.parent / "data" / "SNCA_P37840.fasta"
UNIPROT_ID = "P37840"
PROTEIN_NAME = "SNCA"


def parse_fasta(path: Path) -> str:
    """Read FASTA file and return sequence as a single string."""
    lines = path.read_text().strip().splitlines()
    seq_lines = [l for l in lines if not l.startswith(">")]
    return "".join(seq_lines)


def analyze(sequence: str) -> dict:
    """Return length and amino acid composition (fraction per residue)."""
    length = len(sequence)
    composition = {}
    for aa in sorted(set(sequence)):
        composition[aa] = round(sequence.count(aa) / length, 4)
    return {"length": length, "composition": composition}


def main():
    sequence = parse_fasta(FASTA_PATH)
    result = analyze(sequence)

    print(f"Protein: {PROTEIN_NAME} (UniProt: {UNIPROT_ID})")
    print(f"Length: {result['length']} residues")
    print("\nAmino acid composition:")
    for aa, frac in result["composition"].items():
        count = sequence.count(aa)
        print(f"  {aa}: {count:>3} ({frac:.4f})")

    # Summary for diary.md
    print(f"\n---DIARY---")
    print(f"Protein: {PROTEIN_NAME}")
    print(f"UniProt ID: {UNIPROT_ID}")
    print(f"Length: {result['length']}")
    print(f"Composition: {result['composition']}")


if __name__ == "__main__":
    main()