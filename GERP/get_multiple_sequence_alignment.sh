#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./get_multiple_sequence_alignment.sh <name of reference genome> <input hal file>"
        exit -1
fi



genome=$1
hal=$2


mkdir -p maf/$genome


re='^[0-9]+$'


## Obtain for each chromosome a multiple sequence alignment in MAF format
halStats --bedSequences $genome $hal | cut -f1,3 > maf/$genome/$genome.bed


cat maf/$genome/$genome.bed | while read chromosome length
do
        if [[ $chromosome =~ $re ]] || [[ $chromosome =~ "W" ]] || [[ $chromosome =~ "Z" ]]
        then
                bsub -R'select[mem>180000] rusage[mem=180000]' -M180000 -n 15 -q normal -G rdgroup -J hal2maf -o output_%J -e error_%J hal2maf --refSequence $chromosome \
                --refGenome $genome --noAncestors --onlyOrthologs $hal maf/$genome/$genome.$chromosome.maf
        fi
done
