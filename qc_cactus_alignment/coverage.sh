#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
	echo "Usage: ./coverage.sh <reference.genome> <cactus.alignment>"
	exit -1
fi


REF=$1
HAL=$2


mkdir -p coverage && mkdir -p coverage/$REF

# Obtain a tab delimited genome file that will look like this: <chromName><TAB><chromSize>
halStats --bedSequences $REF $HAL | cut -f1,3 > coverage/$REF/$REF.genome

# Let's now generate 100 random intervals, each 1,000,000 bp long
bedtools random -g coverage/$REF/$REF.genome -l 1000000 -seed 12345 -n 100 > coverage/$REF/$REF.random.intervals

# Let's now get an alignment for each of these 100 random intervals
cat coverage/$REF/$REF.random.intervals | while read contig start end num length strand
do
	bsub -R'select[mem>6000] rusage[mem=6000]' -M6000 -q basement -n 15 -G rdgroup -J coverage -o output_%J -e error_%J hal2maf --refSequence $contig --start $start --length `(expr $end - $start)` --refGenome $REF --onlyOrthologs --noAncestors $HAL coverage/$REF/$REF.random.intervals.$num.maf
done

# Calculate coverage using maf_stream
for num in $(seq 1 100);do /software/team118/maf_stream/target/release/maf_stream coverage $REF coverage/$REF/$REF.random.intervals.$num.maf coverage/$REF/$REF.random.intervals.$num.maf.cov;done

# Plot coverage (boxplot)
# Ideally the species reported in the species_group.tsv file should follow the order of the phylogenetic tree
python3 plot_coverage.py --cov coverage/$REF --refGenome $REF --species species_group.tsv

