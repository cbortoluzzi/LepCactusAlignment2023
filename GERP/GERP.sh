#!/bin/bash



# Author: @cb46




if [ -z $1 ]; then
        echo "Usage: ./GERP.sh <multiple sequence alignment in MAF format> <neutral model>"
        exit -1
fi



MAF=$1
MOD=$2


refGenome=`(basename $MAF .maf) | cut -f1 -d'.'`
chromosome=`(basename $MAF .maf) | cut -f2 -d'.'`


mkdir -p gerp_score/$refGenome


# Filter multiple sequence alignment. We will use default settings: at least 3 species per alignment block and less than 3 species with gapped sequences
python3 filter_alignment.py --maf $MAF --o gerp_score/$refGenome


# Run GERP++ (it has to be installed via conda: https://anaconda.org/bioconda/gerp)
# We will use the --j option to project out the reference genome to avoid any bias
cat $MOD | grep 'TREE' | awk '{print $2}' > nonconserved.$refGenome.mod
gerpcol -t nonconserved.$refGenome.mod -f gerp_score/$refGenome/$refGenome.$chromosome.maf.gerp -e $refGenome -j -x ".4d.rates"


# Assign the GERP score (i.e. rejected substitution score) to each position in the filtered multiple sequence alignment
python3 gerp_score_per_nucleotide_in_alignment.py --maf gerp_score/$refGenome/$refGenome.$chromosome.maf.gerp


