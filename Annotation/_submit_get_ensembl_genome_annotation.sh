#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
  echo "Usage: ./_submit_get_ensembl_genome_annotation.sh <species.list> <species.metadata.json> <output.directory.file>"
  exit -1
fi


SPECIES_LIST=$1
METADATA_JSON=$2
PATHS_FILE=$3


echo "python3 get_ensembl_genome_annotation.py --table $SPECIES_LIST --metadata $METADATA_JSON --path_file $PATHS_FILE"
python3 get_ensembl_genome_annotation.py --table $SPECIES_LIST --metadata $METADATA_JSON --path_file $PATHS_FILE

