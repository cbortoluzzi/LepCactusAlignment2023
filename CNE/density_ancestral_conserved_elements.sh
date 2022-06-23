#!/bin/bash


# Author : @cb46


# Export path
export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares:$PATH


if [ -z $1 ]; then
        echo "Usage: ./density.sh <ancestral conserved elements> <genome file size>"
        exit -1
fi


ANCESTRAL_F=$1
GENOME=$2


mkdir -p ancestral_conserved_elements_density
species_name=$(basename $ANCESTRAL_F | sed 's/ancestral_conserved_elements_//g' | sed 's/.txt//g')

# Calculate density in N-bp window
bedtools makewindows -g $GENOME -w 100000 > ancestral_conserved_elements_density/$species_name.window.bed
# We will consider only ancestral conserved elements that do not contain repeats
cat $ANCESTRAL_F | awk '{if($5 == 0.0)print}' | bedtools coverage -a ancestral_conserved_elements_density/$species_name.window.bed -b stdin > ancestral_conserved_elements_density/$species_name.100kb.bed
rm ancestral_conserved_elements_density/$species_name.window.bed

# Plot density distribution of ancestral conserved elements
species=`echo $species_name | cut -f1,2 -d'_' | sed 's/.*/\u&/'`
/software/R-4.0.3/lib/R/bin/Rscript plot_density_ancestral_conserved_elements.r ancestral_conserved_elements_density/$species_name.100kb.bed $species
