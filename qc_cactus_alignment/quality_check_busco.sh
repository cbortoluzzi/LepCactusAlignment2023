#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./quality_check_busco.sh <name of query species> <input hal file> <name of output directory>"
        exit -1
fi


QUERY=$1
HAL=$2
OUTPUT=$3


# Obtain pairwise comparisons between query and target species
# In our study, the query species was Tinea_trinotella, the most basal of the lepidoptera species
python3 consistency.py --refGenome $QUERY --list_genes final_busco_ids.txt --species_list species_list.tsv --tree supermatrix_datafreeze_080621.treefile.pruned --hal $HAL --o $OUTPUT


# Plot consistent single-copy BUSCO genes per species, as well as inconsistent genes
python3 plot_consistency.py --d busco_quality_check/ --tree supermatrix_datafreeze_080621.treefile.pruned --refGenome $QUERY


