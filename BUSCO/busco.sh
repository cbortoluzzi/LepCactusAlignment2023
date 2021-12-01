#!/bin/bash


# Author: @cb46


# Activate conda environment
conda activate /software/team118/miniconda3/envs/busco5


if [[ -z $1 ]] ; then
  echo "Usage: ./busco.sh <fasta.file> <busco.lineage>"
  exit -1
fi


FASTA=$1
LINEAGE=$2


SPECIES=`echo $FASTA | rev | cut -d"/" -f2 | rev`
DIR=$(dirname $FASTA)

echo "./busco.sh $FASTA $LINEAGE"


echo "busco -i $FASTA -l $LINEAGE -o $SPECIES -c 8 -m genome --offline"
busco -i $FASTA -l $LINEAGE -o $SPECIES -c 8 -m genome --offline


mv $SPECIES vertebrata_odb10_metaeuk
mv vertebrata_odb10_metaeuk $DIR

