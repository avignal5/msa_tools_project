#!/usr/bin/env python3

from msa_tools import read_fasta, save_fasta, read_gff3, save_pairwise_diff_csv, plot_distance_matrix, nucleotide_diversity_fast

#read fasta file
msa = read_fasta("/Users/avignal/Documents/bin/Python/tests/Algeria_Refs_mafft_6_15_Genes_no_ref.fasta")

#distance matrix
names, matrix = msa.distance_matrix()
print(names)
print(matrix)

#plot distance matrix
plot_distance_matrix(names, matrix, "/Users/avignal/Documents/bin/Python/tests/Alg_heatmap.pdf", fontsize=2)
