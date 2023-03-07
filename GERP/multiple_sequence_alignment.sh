#!/bin/bash


# Author: @cb46


if [ -z $1 ]; then
        echo "Usage: ./multiple_sequence_alignment.sh <name of reference genome> <input hal file>"
        exit -1
fi


genome=$1
hal=$2



mkdir -p maf/$genome


echo "Print sequences of given genome in bed format"
halStats --bedSequences $genome $hal | cut -f1,3 | awk '{if($1 < 100 || $1 == "W" || $1 == "Z")print}' | sort -k1,1 -k2,2n > maf/$genome/$genome.txt


echo "Convert hal database to maf"
cat maf/$genome/$genome.txt | while read chromosome length
do
        bsub -R'select[mem>180000] rusage[mem=180000]' -M180000 -n 15 -q basement -G rdgroup -J hal2maf -o output_%J -e error_%J hal2maf --refGenome $genome --refSequence $chromosome --noAncestors --onlyOrthologs $hal maf/$genome/$genome.$chromosome.maf
done


