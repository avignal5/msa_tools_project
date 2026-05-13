#!/usr/bin/env python3

from msa_tools import build_msa_from_getorganelle, read_fasta, save_fasta, MultipleSequenceAlignment

#read fasta file
msa = read_fasta("/Users/avignal/gbigwork/seqapipop/Mitoch_Analysis/Results_Algeria/AAA_rotations/Algeria_Refs_mafft_1.mfa")

#rotate the sequences
msa_r = msa.rotate(6352) #Position of the start of the reference sequence
save_fasta(msa_r, "/Users/avignal/gbigwork/seqapipop/Mitoch_Analysis/Results_Algeria/AAA_rotations/Algeria_Refs_mafft_1_rotate_1.mfa")

