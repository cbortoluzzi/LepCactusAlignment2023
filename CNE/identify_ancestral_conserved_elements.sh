#!/bin/bash


# Author : @cb46


if [ -z $1 ]; then
        echo "Usage: ./identify_ancestral_conserved_elements.sh "
        exit -1
fi



cat genomes_noAncestors.txt



python3 ancestral_conserved_elements_v1.py --maf --refGenome --refSequence --score --length
