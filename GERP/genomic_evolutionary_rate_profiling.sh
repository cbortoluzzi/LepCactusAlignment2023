#!/bin/bash



# Author: @cb46



export PATH=/software/team118/maf_stream/target/release:$PATH



if [ -z $1 ]; then
        echo "Usage: ./genomic_evolutionary_rate_profiling.sh <input maf file>"
        exit -1
fi


maf=$1



## Select sites with at least 3 ungapped sequences
genome=$(basename $maf .maf | cut -f1 -d'.')
mkdir -p bed/$genome
python3 select_sites_with_ungapped_species.py --maf $maf --o bed/$genome
echo "Done with selecting sites"


## Filter input maf file
bed=$(basename $maf .maf)
mkdir -p GERP++/$genome
maf_stream filter --bed bed/$genome/$bed.bed $maf GERP++/$genome/$bed.maf
echo "Done with maf_stream"


## Run GERP++
chromosome=$(basename $maf .maf | cut -f2 -d'.')
re='^[0-9]+$'
if [[ $chromosome =~ $re ]]
then
        cat neutral_mod/$genome/$genome.ave.nonconserved-4d.mod | grep 'TREE:' | sed 's/TREE: //g' > tree.$genome.ave.nh
        gerpcol -t tree.$genome.ave.nh -f GERP++/$genome/$bed.maf -e $genome -j
elif [[ $chromosome =~ "Z" ]]
then
        cat neutral_mod/$genome/$genome.Z.nonconserved-4d.mod | grep 'TREE:' | sed 's/TREE: //g' > tree.$genome.Z.nh
        gerpcol -t tree.$genome.Z.nh -f GERP++/$genome/$bed.maf -e $genome -j
elif [[ $chromosome =~ "W" ]]
then
        cat neutral_mod/$genome/$genome.W.nonconserved-4d.mod | grep 'TREE:' | sed 's/TREE: //g' > tree.$genome.W.nh
        gerpcol -t tree.$genome.W.nh -f GERP++/$genome/$bed.maf -e $genome -j
fi
echo "Done with GERP++"


## Obtain the position for each rejected substitution score in the multiple sequence alignment
python3 gerp_score_per_nucleotide_in_alignment.py --maf GERP++/$genome/$bed.maf


## Filter out sites that do evolve neutrally (GERP score < 0)
cat GERP++/$genome/$bed.maf.rates.bed | awk '{if($6 >= 0)print}' > GERP++/$genome/$bed.maf.rates.conserved.bed


