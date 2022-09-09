#!/bin/bash



# Author: @cb46


mkdir -p blastn

# Generate input in FASTA format
for f in ancestral_conserved_elements/*.bed
do
	genome=$(basename $f .bed)
	cat $f | awk '{print ">"$1":"$2"-"$2+$3"\n"$8}' > blastn/$genome.fa
	# Run Blastn on query
	echo "Running blastn on" $genome
	blastn -db db/outgroups -query blastn/$genome.fa -outfmt 6 -out blastn/$genome.blastn.out
done

