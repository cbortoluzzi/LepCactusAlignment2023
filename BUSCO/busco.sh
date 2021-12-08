#!/bin/bash


# Author: @cb46


# Activate conda environment
conda activate /software/team118/miniconda3/envs/busco5


if [[ -z $1 ]] ; then
  echo "Usage: ./busco.sh <fasta> <busco.lineage>"
  exit -1
fi


FASTA=$1
LINEAGE=$2
DIR=$(dirname $FASTA)


echo "./busco.sh $FASTA $LINEAGE"


echo "busco -i $FASTA -l $LINEAGE -o $DIR/vertebrata_odb10_metaeuk -c 8 -m genome --offline"
busco -i $FASTA -l $LINEAGE -o $DIR/vertebrata_odb10_metaeuk -c 8 -m genome --offline


