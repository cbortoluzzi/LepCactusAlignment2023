#!/bin/bash



# Author: @cb46


export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/gffread:$PATH



if [ -z $1 ]; then
	echo "Usage: ./generate_input.sh <input hal file> <name of reference genome>"
	exit -1
fi



HAL=$1
REF=$2



mkdir -p sequences/$REF && mkdir -p gerp_score/$REF


# Print sequences of a given genome
halStats --bedSequences $REF $HAL | cut -f1,3 > sequences/$REF/$REF.bed


re='^[0-9]+$'
cat sequences/$REF/$REF.bed | while read contig length
do
	if [[ $contig =~ $re ]] || [[ $contig =~ "W" ]] || [[ $contig =~ "Z" ]]
	then
		# Obtain a multiple sequence alignment in MAF format for each contig
		echo "hal2maf --refSequence $contig --refGenome $REF --noAncestors --onlyOrthologs $HAL sequences/$REF/$REF.$contig.maf" > run.hal2maf.$contig.sh
		# Make alignment depth wiggle plot for a genome. By default, this is a count of the number of other unique genomes each base aligns to, including ancestral genomes
		echo "halAlignmentDepth --noAncestors --outWiggle sequences/$REF/$REF.$contig.wig --refSequence $contig $HAL $REF" >> run.hal2maf.$contig.sh
		# Identify consecutive regions with at least 3 ungapped species
		echo "python3 select_sites_with_ungapped_species.py --w sequences/$REF/$REF.$contig.wig --o sequences/$REF" >> run.hal2maf.$contig.sh
		# Filter multiple sequence alignment
		echo "cat sequences/$REF/$REF.$contig.bed | while read chrom begin stop len;do hal2maf --refSequence "'$chrom'" --refGenome $REF --start "'$begin'" --length "'$len'" --noAncestors --onlyOrthologs --append $HAL gerp_score/$REF/$REF.$contig.maf;done" >> run.hal2maf.$contig.sh
	fi
done

chmod +x run.hal2maf.*.sh

