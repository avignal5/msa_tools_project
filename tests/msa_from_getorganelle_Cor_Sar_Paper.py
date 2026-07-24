#!/usr/bin/env python3

from msa_tools import build_msa_from_getorganelle

base_dir = "/Users/avignal/gbigwork/seqapipop/Mitoch_Analysis/Corsica_Sardinia_Paper"
output_mfa = "/Users/avignal/gbigwork/seqapipop/Mitoch_Analysis/Corsica_Sardinia_Paper/test.mfa"
output_csv = "/Users/avignal/gbigwork/seqapipop/Mitoch_Analysis/Corsica_Sardinia_Paper/test.csv"

msa = build_msa_from_getorganelle(base_dir, output_mfa=output_mfa, save_individual=True, reference=None, output_csv=output_csv)

print(msa)
