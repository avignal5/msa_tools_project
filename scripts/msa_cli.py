#!/usr/bin/env python3
"""
msa_cli.py — user-friendly command-line interface for msa_tools

Features:
- Extract regions (single or multiple from BED/list)
- Rotate sequences
- Remove gap-only columns
- Compute pairwise mismatches/gaps
- Export FASTA and CSV statistics
"""

import argparse
import csv
from pathlib import Path
from msa_tools import read_fasta, save_fasta, save_pairwise_diff_csv, MultipleSequenceAlignment, AlignmentError

def parse_bed_file(bed_path: Path):
    """Parse simple BED-like file: chrom/start/end per line"""
    regions = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            start, end = int(parts[1]), int(parts[2])
            regions.append((start + 1, end))  # BED 0-based → 1-based
    return regions


def main():
    parser = argparse.ArgumentParser(description="Friendly MSA CLI")
    parser.add_argument("fasta", type=Path, help="Input aligned FASTA file")
    parser.add_argument("--ref", help="Reference sequence name")
    parser.add_argument("--extract", nargs=2, metavar=("START", "END"), type=int,
                        help="Extract a single region (1-based coordinates)")
    parser.add_argument("--extract-bed", type=Path,
                        help="Extract multiple regions from BED file (0-based)")
    parser.add_argument("--rotate", type=int, help="Rotate sequences by given shift")
    parser.add_argument("--remove-gaps", action="store_true",
                        help="Remove columns containing only gaps")
    parser.add_argument("--pairwise", action="store_true",
                        help="Compute pairwise mismatches and gaps")
    parser.add_argument("--out-fasta", type=Path, help="Save resulting MSA to FASTA")
    parser.add_argument("--out-csv", type=Path, help="Save pairwise stats to CSV")

    args = parser.parse_args()

    try:
        msa = read_fasta(args.fasta, reference=args.ref)

        # Extract regions
        if args.extract:
            start, end = args.extract
            msa = msa.extract_region(start, end)

        if args.extract_bed:
            regions = parse_bed_file(args.extract_bed)
            combined_sequences = {}
            for name in msa.sequences.keys():
                combined_sequences[name] = ""
            for start, end in regions:
                sub_msa = msa.extract_region(start, end)
                for name, seq in sub_msa.sequences.items():
                    combined_sequences[name] += seq
            msa = MultipleSequenceAlignment(combined_sequences, args.ref)

        # Rotate sequences
        if args.rotate:
            msa = msa.rotate(args.rotate)

        # Remove gap-only columns
        if args.remove_gaps:
            msa = msa.remove_gap_columns()

        # Display resulting MSA
        print(f"\nResulting MSA ({len(msa.sequences)} sequences, length={msa.length}):")
        for name, seq in msa.sequences.items():
            print(f">{name}\n{seq}")

        # Pairwise statistics
        stats = None
        if args.pairwise or args.out_csv:
            stats = msa.pairwise_differences()
            if args.pairwise:
                print("\nPairwise mismatches/gaps:")
                for pair, values in stats.items():
                    print(f"{pair[0]} vs {pair[1]}: mismatches={values['mismatches']}, gaps={values['gaps']}")

        # Save outputs
        if args.out_fasta:
            save_fasta(msa, args.out_fasta)
            print(f"\nSaved resulting MSA to {args.out_fasta}")

        if args.out_csv and stats:
            save_pairwise_diff_csv(stats, args.out_csv)
            print(f"Saved pairwise statistics to {args.out_csv}")

    except AlignmentError as e:
        print(f"Error: {e}")
    except FileNotFoundError as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
