#!/bin/bash


# Author : @cb46


if [ -z $1 ]; then
        echo "Usage: ./identify_ancestral_conserved_elements.sh <name of reference sequence within reference genome> <name of reference genome> <input hal file>"
        exit -1
fi



REFSEQ=$1
REFGENOME=$2
HAL=$3


mkdir -p $REFSEQ
cat genomes_noAncestors.txt | while read species
do
        hal2maf --refSequence $REFSEQ --refGenome $REFGENOME --targetGenome $species --onlyOrthologs $HAL $REFSEQ/$species.maf
        # By default, this python script will retain conserved sequences that are at least 50 bp long and are at least 70% identical to the query sequence
        python3 ancestral_conserved_elements_v1.py --maf $REFSEQ/$species.maf --refGenome $REFGENOME --refSequence $REFSEQ
done
