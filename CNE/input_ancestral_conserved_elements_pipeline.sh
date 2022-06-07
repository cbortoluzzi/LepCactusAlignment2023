#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./input_ancestral_conserved_elements_pipeline.sh <name of reference genome> <input hal file>"
        exit -1
fi


REF=$1
HAL=$2



# Print sequences of a given genome in bed format (in our case the genome is the MRCA
halStats --bedSequences $REF $HAL | sort -k3,3 -nr > $REF.sequences.bed
split -l 300 $REF.sequences.bed sequences


# Print the list of genomes in the alignment
halStats --genomes $HAL > genomes.txt
for species in `cat genomes.txt`; do echo $species | grep -v '^Anc';done > genomes_noAncestors.txt
rm genomes.txt
