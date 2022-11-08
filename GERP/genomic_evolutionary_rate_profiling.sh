#!/bin/bash


# Author: @cb46


if [ -z $1 ]; then
        echo "Usage: ./genomic_evolutionary_rate_profiling.sh <name of reference genome>"
        exit -1
fi


genome=$1


mkdir -p bed/$genome && mkdir -p GERP++/$genome


re='^[0-9]+$'
for maf in maf/$genome/*.maf
do
        echo "Select sites with at least N species without gaps"
        python3 select_sites_with_ungapped_species.py --maf $maf --n 3 --o bed/$genome

        echo "Filter multiple sequence alignment"
        maf_stream filter --bed bed/$genome/$(basename $maf .maf).bed $maf GERP++/$genome/$(basename $maf)

        echo "Run GERP++"
        chromosome=$(basename $maf .maf | cut -f2 -d'.')
        if [[ $chromosome =~ $re ]]
        then
                cat nmodel/$genome/$genome.ave.nonconserved-4d.mod | grep 'TREE:' | sed 's/TREE: //g' > nmodel/$genome/$genome.ave.nh
                gerpcol -t nmodel/$genome/$genome.ave.nh -f GERP++/$genome/$(basename $maf) -e $1 -j
                python3 gerp_score_per_nucleotide_in_alignment.py --maf GERP++/$genome/$(basename $maf)
        elif [[ $chromosome == "Z" ]]
        then
                cat nmodel/$genome/$genome.Z.nonconserved-4d.mod | grep 'TREE:' | sed 's/TREE: //g' > nmodel/$genome/$genome.Z.nh
                gerpcol -t nmodel/$genome/$genome.Z.nh -f GERP++/$genome/$(basename $maf) -e $1 -j
                python3 gerp_score_per_nucleotide_in_alignment.py --maf GERP++/$genome/$(basename $maf)
        elif [[ $chromosome == "W" ]]
        then
                cat nmodel/$genome/$genome.W.nonconserved-4d.mod | grep 'TREE:' | sed 's/TREE: //g' > nmodel/$genome/$genome.W.nh
                gerpcol -t nmodel/$genome/$genome.W.nh -f GERP++/$genome/$(basename $maf) -e $1 -j
                python3 gerp_score_per_nucleotide_in_alignment.py --maf GERP++/$genome/$(basename $maf)
        fi
done

rm nmodel/$genome/$genome.*.nh

