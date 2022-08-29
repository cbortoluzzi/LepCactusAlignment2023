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


# Print sequences of given genome in BED format 
halStats --bedSequences $REF $HAL | cut -f1,3 > coverage/$REF/$REF.genome


# Generate 100 random intervals, each 1,000,000 bp long
bedtools random -g coverage/$REF/$REF.genome -l 1000000 -seed 12345 -n 100 > coverage/$REF/$REF.random.intervals


# Convert HAL database to an alignment in multiple alignment format (MAF)
cat coverage/$REF/$REF.random.intervals | while read contig start end num length strand
do
        hal2maf --refSequence $contig --start $start --length 1000000 --refGenome $REF --onlyOrthologs --noAncestors $HAL coverage/$REF/$REF.random.intervals.$num.maf
        # Calculate coverage
        maf_stream coverage $REF coverage/$REF/$REF.random.intervals.$num.maf coverage/$REF/$REF.random.intervals.$num.maf.cov
done

# Plot coverage
python3 plot_coverage.py --species_list species_list.tsv --cov coverage/$REF --t supermatrix_datafreeze_080621.treefile.pruned --o coverage/$REF --refGenome $REF


