#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
	echo "Usage: ./alignment_depth.sh <input hal file>"
	exit -1
fi


hal=$1


mkdir -p coding_sequence_depth


cat species_list.tsv | sed '1d' | while read tol_id class order superfamily family latin_name assembly
do
	 # Path to annotation fasta file
	annotation=/lustre/scratch123/tol/projects/lepidoptera/data/insects/$latin_name/analysis/$tol_id/gene/ensembl/latest/$assembly.ensembl.*.genes.gff3.gz
	# Specify name of species as it appears in the cactus alignment
	species_name=`echo $latin_name | tr '[:upper:]' '[:lower:]'`
	assembly_name=`echo $assembly | sed 's/\./v/g' | sed 's/_//g' | tr '[:upper:]' '[:lower:]'`
	genome=$species_name"_"$assembly_name

	# Copy annotation file
	mkdir -p coding_sequence_depth/$genome
	cp $annotation coding_sequence_depth/$genome/$genome.gff3.gz

	if [ -z "$(ls -A coding_sequence_depth/$genome/)" ]
	then
		rm -r coding_sequence_depth/$genome
	else
		# Obtain unique nonoverlapping coding sequence
		zcat coding_sequence_depth/$genome/$genome.gff3.gz | grep -v '#' | awk '{if($3 == "CDS")print}' | sortBed | bedtools merge -i stdin > coding_sequence_depth/$genome/$genome.merged.CDS.bed
		# Print sequences of given genome in bed format
		halStats --bedSequences $genome $hal | cut -f1,3 | sort -k1,1 -k2,2n | while read chromosome length
		do
			# Run halAlignmentDepth on each protein-coding sequence
			bsub -R'select[mem>160000] rusage[mem=160000]' -M160000 -n 15 -q basement -G rdgroup -J halAlignmentDepth -o output_%J -e error_%J halAlignmentDepth --noAncestors --outWiggle coding_sequence_depth/$genome/$genome.$chromosome.wig --refSequence $chromosome $hal $genome
		done
	fi
done


