#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./identify_ancestral_conserved_elements.sh <reference.genome> <cactus.alignment> <reference.sequence>"
        exit -1
fi


REF=$1
HAL=$2
CONTIG=$3



echo "./identify_aCE.sh $REF $HAL $CONTIG"


mkdir -p $CONTIG
cat genomes_noAncestors.txt | while read targetGenome
do
        echo "Analysing species ..." $targetGenome

        echo "hal2maf --refGenome $REF --refSequence $CONTIG --onlyOrthologs --targetGenomes $targetGenome $HAL $CONTIG/$targetGenome.maf"
        hal2maf --refGenome $REF --refSequence $CONTIG --onlyOrthologs --targetGenomes $targetGenome $HAL $CONTIG/$targetGenome.maf

        echo "python3 ancestral_conserved_elements.py -maf $CONTIG/$targetGenome.maf -refGenome $REF -refSequence $CONTIG"
        python3 ancestral_conserved_elements.py --maf $CONTIG/$targetGenome.maf --refGenome $REF --refSequence $CONTIG
done
