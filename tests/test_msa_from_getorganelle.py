#!/usr/bin/env python3

from msa_tools import build_msa_from_getorganelle

base_dir = "/Users/avignal/gbigwork/seqapipop/Mitoch_Analysis/Results_Sardinia"
output_mfa = "/Users/avignal/gbigwork/seqapipop/Mitoch_Analysis/Results_Sardinia/test.mfa"
output_csv = "/Users/avignal/gbigwork/seqapipop/Mitoch_Analysis/Results_Sardinia/test.csv"

msa = build_msa_from_getorganelle(base_dir, output_mfa=output_mfa, save_individual=False, reference=None, output_csv=output_csv)

print(msa)
