#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./genomes.sh <input hal file>"
        exit -1
fi


HAL=$1


# Print list of species in the alignment
halStats --genomes $HAL | tr -s ' '  '\n' | grep -v '^Anc' | sort > genomes_noAncestors.txt

