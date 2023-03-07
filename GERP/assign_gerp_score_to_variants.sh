#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./assign_gerp_score_to_variants.sh <vcf file> <path to GERP++ output> <name of reference genome>"
        exit -1
fi


vcf=$1
bed=$2
refGenome=$3


# Select heterozygous sites from filtered VCF file
zcat $vcf | grep -v '#' | grep '0/1'| awk 'OFS="\t"{print $1,$2-1,$2,$4,$5}' > $refGenome.heterozygous.sites.bed


# Assign GERP score to heterozygous sites
cat $bed/$refGenome/$refGenome.*.maf.rates.bed | awk '{if($6 >= 5)print}' | bedtools intersect -a stdin -b $refGenome.heterozygous.sites.bed -wa -wb | awk 'OFS="\t"{if($4 == $10 && $2 == $8)print $1,$2,$3,$4,$11,$5,$6}' > $(basename $vcf .gz).GERPscore.bed

rm $refGenome.heterozygous.sites.bed

