#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./quality_check_busco.sh <name of reference genome> <input hal file> <output directory>"
        exit -1
fi


REFGENOME=$1
HAL=$2
OUTPUT=$3


# Obtain pairwise comparisons between query and target species
python3 consistency.py --refGenome $REFGENOME --list_genes final_busco_ids.txt --species_list species_list.tsv --tree supermatrix_datafreeze_080621.treefile.pruned --hal $HAL --o $OUTPUT


# Plot consistent single-copy BUSCO genes per species, as well as inconsistent genes
python3 plot_consistency.py --d busco_quality_check/ --tree supermatrix_datafreeze_080621.treefile.pruned --refGenome $REFGENOME


