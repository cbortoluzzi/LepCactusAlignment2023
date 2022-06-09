#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
                echo "Usage: ./alignment_depth.sh <annotation in GFF (General Feature Format) format> <name of reference genome> <input hal file>"
                exit -1
fi


GFF=$1
REFGENOME=$2
HAL=$3


mkdir -p coding_sequence

# Obtain unique non-overlapping instances of each feature type
python3 feature_selection.py --gff $GFF --feature CDS --o coding_sequence/$REFGENOME


# Count number of unique genomes that aligns to each position in the coding sequence
cat coding_sequence/$REFGENOME/$REFGENOME\_CDS.merged.bed | while read chrom start end feature strand gene
do
        # We set the --plot option to No becasue we don't want to generate a plot for each coding sequence
        # We set the --depth to 1 because we want each position in the coding sequence to map to at least one genome
        python3 alignment_depth.py --hal $HAL --refSequence $chrom --start $start --length `expr $end - $start` --depth 1 --refGenome $REFGENOME --plot No --o coding_sequence/$REFGENOME
done
