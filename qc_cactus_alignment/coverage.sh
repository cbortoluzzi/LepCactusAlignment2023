#!/bin/bash



# Author: @cb46



export PATH=/software/team118/maf_stream/target/release:$PATH



if [ -z $1 ]; then
        echo "Usage: ./coverage.sh <input hal file> <name of genome to calculate coverage on>"
        exit -1
fi



hal=$1
genome=$2


mkdir -p coverage/$genome


# Print sequences of given genome in BED format
halStats --bedSequences $genome $hal | cut -f1,3 | sort -k1,1 -k2,2n | awk '{if($1 < 100 || $1 == "W" || $1 == "Z")print}' > coverage/$genome/$genome.txt


# Generate 100 random intervals, each 1,000,000 bp long
bedtools random -g coverage/$genome/$genome.txt -l 1000000 -seed 12345 -n 100 > coverage/$genome/$genome.random.intervals


# Convert HAL database to a reference-based alignment in multiple alignment format (MAF)
cat coverage/$genome/$genome.random.intervals | while read chromosome start end num length strand
do
	hal2maf --refSequence $chromosome --start $start --length $length --refGenome $genome --onlyOrthologs --noAncestors $hal coverage/$genome/$genome.random.intervals.$num.maf
	# Calculate coverage
	maf_stream coverage $genome coverage/$genome/$genome.random.intervals.$num.maf coverage/$genome/$genome.random.intervals.$num.maf.cov
done


# Plot coverage
python3 plot_coverage.py --t phylogenetic_tree.nh --l 1000000 --c coverage/$genome --f species_list.tsv --refGenome $genome --o coverage/$genome

