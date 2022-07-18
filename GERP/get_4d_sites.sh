#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./get_4d_sites.sh <input maf file>"
        exit -1
fi




maf=$1


chr=$(basename $maf .maf | cut -f2 -d'.')
genome=$(basename $maf .maf | cut -f1 -d'.')

## Obtain 4d sites in SS format
msa_view $maf --in-format MAF --4d --features 4d_sites/$genome/$genome.$chr.gff3 --out-format SS > 4d_sites/$genome/$genome.$chr.codons.ss

## Transform the format of 4-fold generate site form SS to FASTA
msa_view 4d_sites/$genome/$genome.$chr.codons.ss --in-format SS --out-format FASTA --tuple-size 1 > 4d_sites/$genome/$genome.$chr.sites.fa

