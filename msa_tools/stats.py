import numpy as np
from math import comb

def nucleotide_diversity_fast(msa):

    arr = np.array([list(seq) for seq in msa.sequences.values()])
    n = arr.shape[0]

    total = 0.0
    pair_count = 0

    for i in range(n):
        for j in range(i+1, n):

            s1 = arr[i]
            s2 = arr[j]

            mask = (s1 != "-") & (s2 != "-")

            comparable = np.sum(mask)
            if comparable == 0:
                continue

            differences = np.sum(s1[mask] != s2[mask])

            total += differences / comparable
            pair_count += 1

    return total / pair_count if pair_count else 0.0


def sliding_window_pi_fast(msa, window=100, step=10):
    """
    Ultra-fast sliding window nucleotide diversity π.
    Ignores gaps.

    Returns list of (start, end, pi_value)
    """

    arr = np.array([list(seq) for seq in msa.sequences.values()])
    n_seq, L = arr.shape

    # Precompute per-site π
    pi_per_site = np.zeros(L)

    for i in range(L):

        column = arr[:, i]
        non_gap = column[column != "-"]
        n_k = len(non_gap)

        if n_k < 2:
            pi_per_site[i] = 0
            continue

        # count allele frequencies
        values, counts = np.unique(non_gap, return_counts=True)

        total_pairs = comb(n_k, 2)

        diff_pairs = 0
        for c in counts:
            diff_pairs += c * (n_k - c)

        diff_pairs /= 2  # each pair counted twice

        pi_per_site[i] = diff_pairs / total_pairs

    # cumulative sum for fast windows
    cumsum = np.cumsum(pi_per_site)

    results = []

    for start in range(0, L, step):
        end = min(start + window, L)

        if start == 0:
            total = cumsum[end - 1]
        else:
            total = cumsum[end - 1] - cumsum[start - 1]

        length = end - start
        pi_window = total / length if length > 0 else 0

        results.append((start + 1, end, pi_window))

    return results