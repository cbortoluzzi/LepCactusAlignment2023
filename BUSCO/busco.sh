#!/bin/bash


# Author: @cb46


# Activate conda environment
conda activate /software/team118/miniconda3/envs/busco5


if [[ -z $1 ]] ; then
        echo "Usage: ./busco.sh <assembly in fasta format> <busco lineage>"
        exit -1
fi


FASTA=$1
LINEAGE=$2


SPECIES=`echo $FASTA | rev | cut -d'/' -f2 | rev`
OUTPUT=`echo $LINEAGE | rev | cut -d'/' -f2| rev`
DIRECTORY=$(dirname $FASTA)


# Run BUSCO using conda
busco -i $FASTA -l $LINEAGE -o $SPECIES -c 8 -m genome --offline

# Rename and move output
mv $SPECIES $DIRECTORY/$OUTPUT

