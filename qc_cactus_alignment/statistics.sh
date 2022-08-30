#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./statistics.sh <input hal file>"
        exit -1
fi


hal=$1

# Retrieve basic statistics from a HAL database
# Some global information from a HAL file can be quickly obtained using halStats. It will return the number of genomes, their phylogenetic tree, and the size of each array in each genome
halStats $hal > global_stats.txt

# A count of each type of mutation (Insertions, Deletions, Inversions, Duplications, Transpositions, Gap Insertions, Gap Deletions) in each branch of the alignment can be printed out in a table
# The --maxGap option is used to distinguish from small, 'gap' indels and larger indels
# HAL allows gap indels to be nested within larger rearrangements: ex. an inversion with a gap deletion inside would be counted as a single inversion, but an inversion containing a non-gap event would be identified
# as multiple independent inversions. --maxNFraction will prevent rearrangements with missing data as being identified as such.
halSummarizeMutations $hal --maxGap 50 --maxNFraction 0 > summary_mutations_table.txt

