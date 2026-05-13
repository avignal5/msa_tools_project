#!/bin/bash

# extract_genes.bash

/Users/avignal/Documents/bin/Python/msa_tools_project/scripts/msa_cli.py /Users/avignal/gbigwork/seqapipop/Mitoch_Analysis/Results_Algeria/AAA_rotations/Algeria_Refs_mafft_6.mfa \
    --ref NC_001566.1/1-16343 \
    --extract-bed /Users/avignal/Documents/Stats/2025_Mitochondria/Circular_only/mitoch_all.bed \
    --remove-gaps \
    --out-fasta Algeria_Refs_mafft_6_15_Genes.fasta \
    --pairwise \
    --out-csv pairwise.csv
