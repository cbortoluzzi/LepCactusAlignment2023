#!/bin/bash


# Author: @cb46


if [ -z $1 ]; then
        echo "Usage: ./genomic_evolutionary_rate_profiling.sh <name of reference genome>"
        exit -1
fi


genome=$1


mkdir -p GERP++/$genome


cat nmodel/$genome/$genome.ave.nonconserved-4d.mod | grep 'TREE:' | sed 's/TREE: //g' > nmodel/$genome/$genome.ave.nh

for maf in maf/$genome/$genome.*.maf
do
        chromosome=$(basename $maf .maf | cut -f2 -d'.')

        echo "RUN GERP++"
        echo "Retain sites in the alignment with no more than 3 ungapped species"
        python3 select_sites_with_ungapped_species.py --maf $maf --o bed/$genome
        maf_stream filter --bed bed/$genome/$genome.$chromosome.bed $maf GERP++/$genome/$genome.$chromosome.maf

        # Sex chromosomes
        if [[ $chromosome =~ [^0-9] ]]
        then
                echo "Calculate the rejected substitution score for each position in the alignment"
                cat nmodel/$genome/$genome.$chromosome.nonconserved-4d.mod | grep 'TREE:' | sed 's/TREE: //g' > nmodel/$genome/$genome.$chromosome.nh
                # Run gerpcol to calculate rejected substitution score
                gerpcol -t nmodel/$genome/$genome.$chromosome.nh -f GERP++/$genome/$genome.$chromosome.maf -e $genome -j
        # Autosomes
        else
                echo "Calculate the rejected substitution score for each position in the alignment"
                gerpcol -t nmodel/$genome/$genome.ave.nh -f GERP++/$genome/$genome.$chromosome.maf -e $genome -j
        fi
        python3 gerp_score_per_nucleotide_in_alignment.py --maf GERP++/$genome/$genome.$chromosome.maf
done

rm nmodel/$genome/*.nh


