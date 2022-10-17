#!/bin/bash


# Author: @cb46


# Activate conda environment
conda activate /software/team118/miniconda3/envs/busco5


if [[ -z $1 ]] ; then
        echo "Usage: ./busco.sh <assembly in fasta format> <busco lineage> <species name>"
        exit -1
fi


fasta=$1
lineage=$2
species=$3


directory=$(dirname $fasta)


# Run BUSCO using conda
busco -i $fasta -l $lineage -o $species -c 8 -m genome --offline

# Rename and move output
mv $species $directory/$lineage

