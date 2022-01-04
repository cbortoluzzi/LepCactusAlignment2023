#!/bin/bash



# Author: @cb46



echo "Usage: ./concatenate_aCE.sh"


DATE="$(date +'%Y%m%d')"
DIR="conserved_elements_"$DATE


mkdir -p $DIR

cat genomes_noAncestors.txt | while read targetGenome
do
        echo $targetGenome
        conserved_elements=$(echo 'ancestral_conserved_elements_'$targetGenome'.txt')
        # Move ancestral conserved elements ...
        for f in Anc01refChr*/$conserved_elements; do cat $f >> $DIR/$conserved_elements; done
done

