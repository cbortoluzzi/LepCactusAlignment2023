#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
	echo "Usage: ./get_unique_cds.sh <name of reference genome>"
	exit -1
fi


genome=$1


mkdir -p 4d_sites/$genome


# Merge overlapping coding sequences
cat maf/$genome/$genome.bed | while read chromosome length
do
	zcat gene/$genome/$genome.gff3.gz | grep -v '#' | awk 'OFS="\t"{if($3 == "CDS" && $1 == "'$chromosome'")print}' | sortBed | bedtools merge -i stdin -d -1 | awk 'OFS="\t"{print "'$genome'."$1,"ensembl","CDS",$2,$3,".","+",".","ID=CDS"}' > 4d_sites/$genome/$genome.$chromosome.gff3
done

