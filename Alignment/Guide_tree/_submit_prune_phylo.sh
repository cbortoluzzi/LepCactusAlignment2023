#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./_submit_prune_phylo.sh <phylogenetic.tree> <species.list> <pruned.phylogenetic.tree>"
        exit -1
fi


TREE=$1
SPECIES_LIST=$2
PRUNED_PHYLO=$3


# Prune phylogenetic tree
echo "python3 prune_phylo.py --tree $TREE --table $SPECIES_LIST --pruned $PRUNED_PHYLO"
python3 prune_phylo.py --tree $TREE --table $SPECIES_LIST --pruned $PRUNED_PHYLO

# Assign color to species based on group. This file can be used in iTOL: https://itol.embl.de/upload.cgi
echo "python3 annotate_phylo.py --table $SPECIES_LIST"
python3 annotate_phylo.py --table $SPECIES_LIST


