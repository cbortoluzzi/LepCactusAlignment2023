#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./input_ancestral_conserved_elements_pipeline.sh <name of reference genome> <input hal file>"
        exit -1
fi


REF=$1
HAL=$2



# Print sequences of a given genome in bed format (in our case the genome is the MRCA). We will retain only sequences with a minimum length of 50 bp
halStats --bedSequences $REF $HAL | awk '{if($3 >= 50)print}' | sort -k3,3 -nr > $REF.sequences.50bp.bed
split -l 300 $REF.sequences.50bp.bed sequences


# Print the list of genomes in the alignment
halStats --genomes $HAL > genomes.txt
for species in `cat genomes.txt`; do echo $species | grep -v '^Anc';done > genomes_noAncestors.txt
rm genomes.txt
