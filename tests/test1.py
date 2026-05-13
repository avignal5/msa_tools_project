#!/usr/bin/env python3

from msa_tools import read_fasta, save_fasta, read_gff3, save_pairwise_diff_csv, plot_distance_matrix, nucleotide_diversity_fast

#read fasta file
msa = read_fasta("/Users/avignal/Documents/bin/Python/msa_tools_project/examples/test2.fasta", reference="seq1")
print(msa)

#save the file
save_fasta(msa, "/Users/avignal/Documents/bin/Python/msa_tools_project/examples/out.fasta")

#rotate the sequences
msa_r = msa.rotate(3)
save_fasta(msa_r, "/Users/avignal/Documents/bin/Python/msa_tools_project/examples/out_rotated.fasta")

#compute pairwise differences
pairwise_differences = msa.pairwise_differences()
print(pairwise_differences)
save_pairwise_diff_csv(pairwise_differences, "/Users/avignal/Documents/bin/Python/msa_tools_project/examples/pairwise.diff")

#distance matrix
names, matrix = msa.distance_matrix()
print(names)
print(matrix)

#plot distance matrix
plot_distance_matrix(names, matrix, "/Users/avignal/Documents/bin/Python/msa_tools_project/examples/heat.png")

#consensus sequence
consensus_seq = msa.consensus()
print(consensus_seq)

#nucleotide diversity pi
print(msa.nucleotide_diversity())
print(nucleotide_diversity_fast(msa))

#Remove columns with gaps in all sequences
msa_no_gaps = msa.remove_gap_columns()
save_fasta(msa_no_gaps, "/Users/avignal/Documents/bin/Python/msa_tools_project/examples/out_no_gaps.fasta")

# Extraction BED
#regions = read_bed("/Users/avignal/Documents/bin/Python/msa_tools_project/examples/regions.bed")
#for start, end in regions:
#    sub_msa = msa.extract_region(start, end)
#    print(sub_msa)

# Extraction GFF3
#genes = read_gff3("/Users/avignal/Documents/bin/Python/msa_tools_project/examples/annotations.gff3")
#gene_msa = msa.extract_features(genes)
#print(gene_msa)
