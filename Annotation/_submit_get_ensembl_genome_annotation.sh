#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
  echo "Usage: ./_submit_get_ensembl_genome_annotation.sh <species.list> <species.metadata.json> <output.directory.file>"
  exit -1
fi


SPECIES=$1
METADATA=$2
PATH=$3


echo "python3 get_ensembl_genome_annotation.py --table $SPECIES --metadata $METADATA --path_file $PATH"
python3 get_ensembl_genome_annotation.py --table $SPECIES --metadata $METADATA --path_file $PATH

