#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./phyloFit.sh <name of reference genome>"
        exit -1
fi



genome=$1


mkdir -p neutral_mod/$genome

for maf in maf/$genome/*.maf;
do
        chr=$(basename $maf .maf | cut -f2 -d'.')
        # Estimate the nonconserved (neutral) model for each chromosome
        phyloFit --tree `cat tree_topology.nh` --subst-mod REV --EM --msa-format FASTA 4d_sites/$genome/$genome.$chr.sites.fa --out-root 4d_sites/$genome/$genome.$chr.nonconserved-4d
done

## Obtain the average nonconserved model for the autosomes
re='^[0-9]+$'
autosomes=`cat maf/$genome/$genome.bed | while read contig length;do if [[ $contig =~ $re ]]; then echo $contig; fi;done`
mod=`for chr in $autosomes;do ls 4d_sites/$genome/$genome.$chr.nonconserved-4d.mod;done`
phyloBoot --read-mods $mod --output-average neutral_mod/$genome/$genome.ave.nonconserved-4d.mod


## Move nonconserved model of sex chromosomes to the right folder
mv 4d_sites/$genome/$genome.Z.nonconserved-4d.mod neutral_mod/$genome
mv 4d_sites/$genome/$genome.W.nonconserved-4d.mod neutral_mod/$genome

