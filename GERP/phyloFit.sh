#!/bin/bash



# Author: @cb46


export PATH=/software/team118/phast/bin:$PATH



if [ -z $1 ]; then
	echo "Usage: ./phyloFit.sh <input hal file> <name of reference genome>"
	exit -1
fi



HAL=$1
REF=$2



mkdir -p msa_view/$REF && mkdir -p neutral_model/$REF


# Extract 4d sites: since we have an alignment for each contig, we will estimate a nonconserved (neutral) model for each of them and then calculate the average nonconserved model
# A separate nonconserved model will be estimated for the W and Z chromosome
for f in sequences/$REF/*.maf
do
	chr=$(basename $f .maf | cut -f2 -d'.')
	# Obtain 4d sites in SS format
	msa_view sequences/$REF/$REF.$chr.maf --in-format MAF --4d --features annotation/$REF/$REF.$chr.CDS.gff3 --out-format SS > msa_view/$REF/$REF.$chr.codons.ss
	msa_view msa_view/$REF/$REF.$chr.codons.ss --in-format SS --out-format SS --tuple-size 1 > msa_view/$REF/$REF.$chr.sites.ss
	# Estimate the nonconserved (neutral) model for each chromosome
	phyloFit --tree `cat tree_topology.nh` --subst-mod REV --EM --msa-format SS msa_view/$REF/$REF.$chr.sites.ss --out-root msa_view/$REF/$REF.$chr.nonconserved-4d
done


# Obtain the average nonconserved model
re='^[0-9]+$'
autosomes=`cat sequences/$REF/$REF.bed | while read contig size;do if [[ $contig =~ $re ]]; then echo $contig; fi;done`
mod=`for chr in $autosomes;do ls msa_view/$REF/$REF.$chr.nonconserved-4d.mod;done`
phyloBoot --read-mods $mod --output-average neutral_model/$REF/$REF.ave.nonconserved-4d.mod

mv msa_view/$REF/$REF.W.nonconserved-4d.mod neutral_model/$REF
mv msa_view/$REF/$REF.Z.nonconserved-4d.mod neutral_model/$REF

rm -r msa_view
