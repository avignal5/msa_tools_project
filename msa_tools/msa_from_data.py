#!/usr/bin/env python3

from pathlib import Path
from msa_tools import MultipleSequenceAlignment, save_fasta
import csv

def build_msa_from_getorganelle(base_dir, output_mfa=None, save_individual=False, reference=None, output_csv=None):
    """
    Build a MultipleSequenceAlignment from sequences in getorganelle subfolders.

    Parameters
    ----------
    base_dir : str or Path
        Base directory containing subfolders with sample_name/sequence.fa
        Expected structure: base_dir/seq_name/sequence.fa
    output_mfa : str or Path, optional
        Path to write combined MSA (.mfa)
    save_individual : bool
        If True, save renamed seq_name.fa in each folder
    reference : str, optional
        Name of the reference sequence for MSA object

    Returns
    -------
    msa : MultipleSequenceAlignment

    Usage example
    -------
    base_dir = "/data"
    output_mfa = "/data/msa.mfa"
    msa = build_msa_from_getorganelle(base_dir, output_mfa=output_mfa, save_individual=True, reference=None)
    print(msa)
    """

    base_dir = Path(base_dir)
    sequences = {}

    data = {}
    for fa_path in sorted(base_dir.glob("*/*.fasta")):

        seq_name = fa_path.parent.name  # folder name

        with open(fa_path) as f:
            seq_lines = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    old_name = line[1:].strip()
                    continue
                seq_lines.append(line)

        sequence = "".join(seq_lines)
        sequences[seq_name] = sequence

        data[seq_name] = (len(sequence), old_name)

        # Save individual renamed file
        if save_individual:
            renamed_path = fa_path.parent / f"{seq_name}.fa"
            with open(renamed_path, "w") as rf:
                rf.write(f">{seq_name}\n{sequence}\n")

    # Build MSA object : not possible, as sequences have different lengths
    #msa = MultipleSequenceAlignment(sequences, reference=reference)

    # Write combined MFA if requested
    if output_mfa:
        output_mfa = Path(output_mfa)
        with open(output_mfa, "w") as out_f:
            for name, seq in sequences.items():
                out_f.write(f">{name}\n{seq}\n")
        print(f"Combined MSA saved to {output_mfa}")

    # Write data on the sequence lengths and assembly paths
    if output_csv:
        with open(output_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Sample", "Length", "Path"])  # header
            for key, (length, old_name) in data.items():
                writer.writerow([key, length, old_name])


#    return msa