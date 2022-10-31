#!/bin/bash


# Author: @cb46


export PATH=/software/team118/phast/bin:$PATH


if [ -z $1 ]; then
        echo "Usage: ./neutral_model.sh <name of reference genome>"
        exit -1
fi


genome=$1


echo "./neutral_model.sh $genome"


mkdir -p 4d_sites/$genome && mkdir -p nmodel/$genome


for bed in maf/$genome/*.bed
do
        echo "Merge overlapping protein coding sequences"
        zcat gene/$genome/$genome.gff3.gz | grep -v '#' | awk '{if($3 == "CDS")print}' | bedtools intersect -b $bed -a stdin | sortBed | mergeBed | \
        awk 'OFS="\t"{print "'$1'."$1,"ensembl","CDS",$2,$3,".","+",".","ID=CDS"}' > 4d_sites/$genome/$genome.$(basename $bed .bed).uniqueCDS.gff3

        echo "Estimate neutral (or nonconserved model)"
        msa_view maf/$genome/$genome.$(basename $bed .bed).maf --in-format MAF --4d --features 4d_sites/$genome/$genome.$(basename $bed .bed).uniqueCDS.gff3 --out-format SS > $output.codons.ss
        msa_view 4d_sites/$genome/$(basename $maf .maf).codons.ss --in-format SS --out-format FASTA --tuple-size 1 > 4d_sites/$genome/$(basename $maf .maf).sites.fa
        phyloFit --tree `cat tree_topology.nh` --subst-mod REV --EM --msa-format FASTA 4d_sites/$genome/$(basename $maf .maf).sites.fa --out-root 4d_sites/$genome/$(basename $maf .maf).nonconserved-4d
done


# Average nonconserved model for the autosomes
re='^[0-9]+$'
autosomes=`cat maf/$genome/$genome.txt | while read contig length;do if [[ $contig =~ $re ]]; then echo $contig; fi;done`
mod=`for chr in $autosomes;do ls 4d_sites/$genome/$genome.$chr.*.nonconserved-4d.mod;done`
phyloBoot --read-mods $mod --output-average nmodel/$genome/$genome.ave.nonconserved-4d.mod

mod=`ls 4d_sites/$genome/$genome.Z.*.nonconserved-4d.mod`
phyloBoot --read-mods $mod --output-average nmodel/$genome/$genome.Z.nonconserved-4d.mod

mod=`ls 4d_sites/$genome/$genome.W.*.nonconserved-4d.mod`
phyloBoot --read-mods $mod --output-average nmodel/$genome/$genome.W.nonconserved-4d.mod

