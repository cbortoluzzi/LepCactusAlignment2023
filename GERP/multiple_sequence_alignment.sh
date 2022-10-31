#!/bin/bash


# Author: @cb46


if [ -z $1 ]; then
        echo "Usage: ./multiple_sequence_alignment.sh <name of reference genome> <input hal file>"
        exit -1
fi


genome=$1
hal=$2


echo "./multiple_sequence_alignment.sh $genome $hal"


mkdir -p maf/$genome


re='^[0-9]+$'

echo "Print sequences of given genome in bed format"
halStats --bedSequences $genome $hal | cut -f1,3 | sort -k1,1 -k2,2n | while read chromosome length
do
        if [[ $chromosome =~ $re ]] || [[ $chromosome =~ "W" ]] || [[ $chromosome =~ "Z" ]]
        then
                printf "%s\t%s\n" $chromosome $length >> maf/$genome/$genome.txt
        fi
done


echo "Makes adjacent or sliding windows across a genome or BED file"
bedtools makewindows -g maf/$genome/$genome.txt -w 2000000 | while read chromosome start end;do printf "%s\t%s\t%s\n" $chromosome $start $end > maf/$genome/$chromosome.$start.$end.bed;done


echo "Convert hal database to maf"
for bed in maf/$genome/*.bed
do
        bsub -R'select[mem>180000] rusage[mem=180000]' -M180000 -n 15 -q basement -G rdgroup -J hal2maf -o output_%J -e error_%J hal2maf --refGenome $genome \
        --noAncestors --onlyOrthologs --refTargets $bed $hal maf/$genome/$genome.$(basename $bed .bed).maf
done

