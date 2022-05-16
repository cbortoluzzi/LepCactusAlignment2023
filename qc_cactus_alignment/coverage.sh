#!/bin/bash



# Author: @cb46


export PATH=/software/team118/maf_stream/target/release:$PATH


if [ -z $1 ]; then
        echo "Usage: ./coverage.sh <input hal file> <genome to calculate coverage on>"
        exit -1
fi


HAL=$1
REF=$2


mkdir -p coverage/$REF


# Obtain a tab delimited genome file that will look like this: <chromName><TAB><chromSize>
halStats --bedSequences $REF $HAL | cut -f1,3 > coverage/$REF/$REF.genome


# Let's now generate 100 random intervals, each 1,000,000 bp long
bedtools random -g coverage/$REF/$REF.genome -l 1000000 -seed 12345 -n 100 > coverage/$REF/$REF.random.intervals


# Let's now get an alignment for each of these 1,000 random intervals
cat coverage/$REF/$REF.random.intervals | while read contig start end num length strand
do
        echo "Generating an alignment for random interval number" $num
        hal2maf --refSequence $contig --start $start --length `(expr $end - $start)` --refGenome $REF --onlyOrthologs --noAncestors $HAL coverage/$REF/$REF.random.intervals.$num.maf
        # Calculate coverage
        maf_stream coverage $REF coverage/$REF/$REF.random.intervals.$num.maf coverage/$REF/$REF.random.intervals.$num.maf.cov
done

# Plot coverage
python3 plot_coverage.py --species_list species_list_012022_v2.tsv --cov coverage/$REF --t supermatrix_datafreeze_080621.treefile.pruned --o coverage/$REF --refGenome $REF

