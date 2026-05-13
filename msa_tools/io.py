from pathlib import Path
from .alignment import MultipleSequenceAlignment
import csv

def read_fasta(filepath: Path, reference: str = None):

    filepath = Path(filepath)

    sequences = {}
    name = None
    seq_chunks = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(seq_chunks)
                name = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if name:
            sequences[name] = "".join(seq_chunks)

    return MultipleSequenceAlignment(sequences, reference)

def save_fasta(msa: MultipleSequenceAlignment, output_path: Path):

    output_path = Path(output_path)

    with open(output_path, "w") as f:
        for name, seq in msa.sequences.items():
            f.write(f">{name}\n{seq}\n")

def save_pairwise_diff_csv(stats: dict, output_path: Path):

    """Save pairwise differences and gaps as CSV"""

    output_path = Path(output_path)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Seq1", "Seq2", "Mismatches", "Gaps", "identity"])
        for (s1, s2), values in stats.items():
            writer.writerow([s1, s2, values["mismatches"], values["gaps"], values["identity"]])

def read_bed(filepath: Path):
    """
    Parse BED file and return list of (start, end) tuples (1-based)
    """

    filepath = Path(filepath)

    regions = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            start, end = int(parts[1]), int(parts[2])
            regions.append((start+1, end))  # convert 0-based BED → 1-based
    return regions