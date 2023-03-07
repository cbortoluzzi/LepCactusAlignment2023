#!/bin/bash


# Author: @cb46


export PATH=/software/team118/phast/bin:$PATH


if [ -z $1 ]; then
        echo "Usage: ./neutral_model.sh <name of reference genome>"
        exit -1
fi


genome=$1


mkdir -p 4d_sites/$genome && mkdir -p nmodel/$genome


cat maf/$genome/$genome.txt | cut -f1 | while read chromosome
do
        echo "Analysing chromosome" $chromosome

        echo "Obtain protein-coding sequences"
        # The descript assumes that there is already a folder for each species, called geneset, with an annotation file in GFF3 format
        zcat geneset/$genome/$genome.gff3.gz | grep -v '#' | awk 'OFS="\t"{if($3 == "CDS" && $1 == "'$chromosome'")print}' | sortBed | mergeBed -i stdin | \
        awk 'OFS="\t"{print "'$genome.'""'$chromosome'","ensembl","CDS",$2,$3,".","+",".","ID=CDS"}' > geneset/$genome/$genome.$chromosome.gff3

        echo "Estimate neutral (or nonconserved) model"
        msa_view maf/$genome/$genome.$chromosome.maf --in-format MAF --4d --features geneset/$genome/$genome.$chromosome.gff3 --out-format SS > 4d_sites/$genome/$genome.$chromosome.codons.ss
        msa_view 4d_sites/$genome/$genome.$chromosome.codons.ss --in-format SS --out-format FASTA --tuple-size 1 > 4d_sites/$genome/$genome.$chromosome.sites.fa
        phyloFit --tree `cat tree_topology.nh` --subst-mod REV --EM --msa-format FASTA 4d_sites/$genome/$genome.$chromosome.sites.fa --out-root 4d_sites/$genome/$genome.$chromosome.nonconserved-4d
done


echo "Average neutral (or nonconserved) model for autosomes"
autosomes=`cat maf/$genome/$genome.txt | awk '{if($1 < 100)print $1}' | sort -k1,1n`
mod=`for chr in $autosomes;do ls 4d_sites/$genome/$genome.$chr.nonconserved-4d.mod;done`
phyloBoot --read-mods $mod --output-average nmodel/$genome/$genome.ave.nonconserved-4d.mod
mv 4d_sites/$genome/$genome.X.nonconserved-4d.mod nmodel/$genome
mv 4d_sites/$genome/$genome.Y.nonconserved-4d.mod nmodel/$genome


