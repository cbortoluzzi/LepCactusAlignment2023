#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./identify_ancestral_conserved_elements.sh <name of reference genome> <input hal file> <name of reference sequence within reference genome>"
        exit -1
fi


REF=$1
HAL=$2
CONTIG=$3


mkdir -p $CONTIG
targetGenomes=`awk '{print $1}' genomes_noAncestors.txt | paste -s -d, -`

# Get alignment
hal2maf --refGenome $REF --refSequence $CONTIG --onlyOrthologs --targetGenomes $targetGenomes $HAL $CONTIG/$CONTIG.maf
# Identify ancestral conserved elements
python3 ancestral_conserved_elements_v2.py --maf $CONTIG/$CONTIG.maf --refGenome $REF --refSequence $CONTIG --score 70 --length 50

