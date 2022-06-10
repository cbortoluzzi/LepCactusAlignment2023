#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./quality_check_busco.sh <name of species to use as query> <input hal file> <name of output directory>"
        exit -1
fi


QUERY=$1
HAL=$2
OUTPUT=$3


# Obtain pairwise comparisons between query (in our case Tinea trinotella) and target species
python3 consistency.py --refGenome $QUERY --list_genes final_busco_ids.txt --species_list species_list.tsv --tree supermatrix_datafreeze_080621.treefile.pruned --hal $HAL --o $OUTPUT


# Plot consistent of single-copy BUSCO genes per species, as well as inconsistent genes
python3 plot_consistency.py --d busco_quality_check/ --tree supermatrix_datafreeze_080621.treefile.pruned --refGenome $QUERY --species_list species_list.tsv --o busco_quality_check


