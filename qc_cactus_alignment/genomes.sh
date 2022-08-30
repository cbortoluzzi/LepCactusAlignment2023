#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./genomes.sh <input hal file>"
        exit -1
fi


hal=$1


# Print list of species in the alignment
halStats --genomes $hal | tr -s ' '  '\n' | grep -v '^Anc' | sort > genomes_noAncestors.txt

