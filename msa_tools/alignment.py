from typing import Dict, Tuple
from collections import Counter
import numpy as np
from .exceptions import AlignmentError

class MultipleSequenceAlignment:
    """
    Representation of a multiple sequence alignment (MSA).
    """

    def __init__(self, sequences: Dict[str, str], reference: str = None):

        if not sequences:
            raise AlignmentError("No sequences provided.")

        lengths = {len(seq) for seq in sequences.values()}
        if len(lengths) != 1:
            raise AlignmentError("All sequences must have same length.")

        self.sequences = sequences
        self.length = lengths.pop()
        self.reference = reference

        if reference and reference not in sequences:
            raise AlignmentError(f"Reference '{reference}' not found.")

    def __repr__(self):
        return f"<MSA: {len(self.sequences)} sequences, length={self.length}>"

    # ---------------------------------------------------------
    # Coordinate mapping
    # ---------------------------------------------------------

    def _build_coordinate_map(self):       # the underscore : internal function
        if not self.reference:
            raise AlignmentError("Reference required.")

        ref_seq = self.sequences[self.reference]
        mapping = {}
        ungapped = 0

        for aln_pos, base in enumerate(ref_seq):
            if base != "-":
                ungapped += 1
                mapping[ungapped] = aln_pos

        return mapping

    # ---------------------------------------------------------
    # A subset of sequences
    # ---------------------------------------------------------

    def subset(self, names):
        """Return a new alignment containing only the selected sequences."""

        missing = set(names) - set(self.sequences)
        if missing:
            raise AlignmentError(
                f"Unknown sequences: {', '.join(sorted(missing))}"
            )

        reference = self.reference if self.reference in names else None

        return MultipleSequenceAlignment(
            {n: self.sequences[n] for n in names},
            reference=reference,
        )

    # ---------------------------------------------------------
    # Reverse complement all or part of a sequence
    # ---------------------------------------------------------

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        table = str.maketrans(
            "ACGTNacgtn-",
            "TGCANtgcan-",
        )
        return seq.translate(table)[::-1]

    def reverse_complement(
        self,
        sequence_name: str,
        start: int = None,
        end: int = None,
        coordinates="alignment",
    ):
        """
        Reverse-complement an entire sequence or a region and return a new MSA.

        Parameters
        ----------
        sequence_name : str
            Sequence to modify.
        start, end : int, optional
            Reference coordinates of the region to reverse-complement.
            If omitted, the entire sequence is reverse-complemented.
        """

        if sequence_name not in self.sequences:
            raise AlignmentError(
                f"Unknown sequence '{sequence_name}'."
            )

        new_sequences = self.sequences.copy()
        seq = new_sequences[sequence_name]

        if start is None and end is None:
            new_sequences[sequence_name] = self._reverse_complement(seq)

        else:

            if coordinates == "reference":
                mapping = self._build_coordinate_map()

                if start not in mapping or end not in mapping:
                    raise AlignmentError(
                        "Coordinates outside reference."
                    )

                start = mapping[start]
                end = mapping[end] + 1

            elif coordinates == "alignment":
                start = start
                end = end + 1

            else:
                raise AlignmentError(
                    "coordinates must be 'alignment' or 'reference'"
                )

            fragment = seq[start:end]

            new_sequences[sequence_name] = (
                seq[:start]
                + self._reverse_complement(fragment)
                + seq[end:]
            )

        return MultipleSequenceAlignment(
            new_sequences,
            self.reference
        )



    # ---------------------------------------------------------
    # Extraction of specific regions
    # ---------------------------------------------------------

    def extract_region(self, start: int, end: int):

        mapping = self._build_coordinate_map()

        if start not in mapping or end not in mapping:
            raise AlignmentError("Coordinates outside reference.")

        aln_start = mapping[start]
        aln_end = mapping[end] + 1

        new_sequences = {
            name: seq[aln_start:aln_end]
            for name, seq in self.sequences.items()
        }

        return MultipleSequenceAlignment(new_sequences, self.reference)

    # ---------------------------------------------------------
    # Pairwise
    # ---------------------------------------------------------

    def pairwise_differences(self):

        """
        Pairwise difference between all sequences in the alignment.
        Counts mismatches and gaps separately.
        """

        names = list(self.sequences.keys())
        results = {}

        for i in range(len(names)):
            for j in range(i + 1, len(names)):

                s1 = self.sequences[names[i]]
                s2 = self.sequences[names[j]]

                mismatches = 0
                gaps = 0

                for a, b in zip(s1, s2):
                    if a == "-" or b == "-":
                        if a != b:
                            gaps += 1
                    elif a != b:
                        mismatches += 1
                matches = sum(1 for a,b in zip(s1,s2) if a==b)
                length = len(s1)
                identity = matches / length * 100
                results[(names[i], names[j])] = {
                    "mismatches": mismatches,
                    "gaps": gaps,
                    "identity": identity
                }

        return results

    def distance_matrix(self):
        """
        Return NxN numpy array of mismatch counts
        """
        names = list(self.sequences.keys())
        N = len(names)
        matrix = np.zeros((N, N), dtype=int)
        for i in range(N):
            for j in range(i+1, N):
                s1, s2 = self.sequences[names[i]], self.sequences[names[j]]
                mismatches = sum(1 for a,b in zip(s1,s2) if a!=b)
                matrix[i,j] = mismatches
                matrix[j,i] = mismatches
        return names, matrix

    # ---------------------------------------------------------
    # Consensus sequence
    # ---------------------------------------------------------

    def consensus(self):
        """
        Return consensus sequence of MSA (most frequent base per column)
        """
        seqs = list(self.sequences.values())
        length = self.length
        consensus_seq = ""
        for i in range(length):
            bases = [s[i] for s in seqs]
            c = Counter(bases)
            consensus_seq += c.most_common(1)[0][0]
        return consensus_seq

    # ---------------------------------------------------------
    # Nucleotide diversity (pi): CHECK THE FORMULA: does not seem correct here and avoid counting gaps as differences SEE SCRIPTS IN stats.py
    # ---------------------------------------------------------

    def nucleotide_diversity(self):
        """
        Calculate π (average pairwise nucleotide diversity)
        """
        names = list(self.sequences.keys())
        N = len(names)
        total = 0
        count = 0
        for i in range(N):
            for j in range(i+1, N):
                s1, s2 = self.sequences[names[i]], self.sequences[names[j]]
                total += sum(1 for a,b in zip(s1,s2) if a!=b)
                count += len(s1)
        return total / count

    """
    # ---------------------------------------------------------
    # For stats on sliding windows HAS NOT BEEN EVALUATED. SEE SCRIPTS IN stats.py
    # ---------------------------------------------------------

    def sliding_window(self, window=50, step=10, metric="diversity"):
    """
        #Apply metric over sliding windows
        #metric: "diversity" (π) or "identity"
    """
        results = []
        L = self.length
        for start in range(0, L, step):
            end = min(start + window, L)
            sub_sequences = {n: s[start:end] for n,s in self.sequences.items()}
            sub_self = type(self)(sub_sequences, self.reference)
            if metric=="diversity":
                val = nucleotide_diversity(sub_self)
            elif metric=="identity":
                val = np.mean(list(pairwise_identity(sub_self).values()))
            else:
                val = None
            results.append((start+1, end, val))
        return results
    """

    # ---------------------------------------------------------
    # Rotation
    # ---------------------------------------------------------

    def rotate(self, shift: int):
        
        """Rotate all sequences in the alignment by `shift` columns and return a new MSA."""

        shift = shift % self.length

        new_sequences = {
            name: seq[shift:] + seq[:shift]
            for name, seq in self.sequences.items()
        }

        return MultipleSequenceAlignment(new_sequences, self.reference)

    # ---------------------------------------------------------
    # Remove gap-only columns
    # ---------------------------------------------------------

    def remove_gap_columns(self):
        
        """Removes columns for which there are gaps in all sequences"""

        array = np.array([list(seq) for seq in self.sequences.values()])
        keep = ~(np.all(array == "-", axis=0))

        new_array = array[:, keep]

        new_sequences = {
            name: "".join(new_array[i])
            for i, name in enumerate(self.sequences.keys())
        }

        return MultipleSequenceAlignment(new_sequences, self.reference)

    # ---------------------------------------------------------
    # Remove gap-only columns
    # ---------------------------------------------------------


    def extract_features(self, features):
        """
        Extract and concatenate regions corresponding to a list of GFFFeature objects.

        Parameters
        ----------
        features : list of GFFFeature

        Returns
        -------
        MultipleSequenceAlignment
        """
        if not self.reference:
            raise AlignmentError("Reference required for GFF extraction.")

        combined_sequences = {name: "" for name in self.sequences.keys()}

        for feature in features:
            sub_msa = self.extract_region(feature.start, feature.end)

            for name, seq in sub_msa.sequences.items():
                combined_sequences[name] += seq

        return MultipleSequenceAlignment(combined_sequences, self.reference)