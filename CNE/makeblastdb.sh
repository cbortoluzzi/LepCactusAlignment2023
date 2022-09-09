#!/bin/bash



# Author: @cb46



# Generate multiple sequence alignment in FASTA format
python3 mfa.py --i list_outgroups.tsv --p /lustre/scratch123/tol/projects/darwin/data/insects/ --o db

# Generate database
makeblastdb -in db/outgroups.fasta -out db/outgroups -parse_seqids -dbtype nucl -input_type fasta

