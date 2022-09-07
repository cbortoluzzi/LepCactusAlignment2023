#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./input_ancestral_conserved_elements_pipeline.sh <name of reference genome> <input hal file>"
        exit -1
fi


ref=$1
hal=$2



# Print sequences of a given genome in bed format (in our case the genome is the MRCA
halStats --bedSequences $ref $hal | sort -k3,3 -nr > $ref.sequences.bed
split -l 300 $ref.sequences.bed sequences_


# Print the list of genomes in the alignment
halStats --genomes $hal > genomes.txt
for species in `cat genomes.txt`; do echo $species | grep -v '^Anc';done > genomes_noAncestors.txt
rm genomes.txt

