#!/bin/bash



# Author: @cb46


if [ -z $1 ]; then
	echo "Usage: ./qc_busco.sh <species.list> <cactus.alignment>"
	exit -1
fi


FILE=$1
HAL=$2

# The input file should contain information on - in order - assembly, tol_id, species latin name, and species name as it appears in the cactus alignment

mkdir -p busco

# Count number of genomes
num_genomes=`cat $FILE | wc -l`

# Identify complete, single copy busco genes that are shared among all species in the alignment
cat $FILE | while read assembly tol_id species genome;do cat /lustre/scratch123/tol/projects/lepidoptera/data/insects/$species/analysis/$tol_id/busco/mm49_lepidoptera_odb10_metaeuk/run_lepidoptera_odb10/full_table.tsv | grep -v '^#' | awk '$2=="Complete" {print $1}'; done >> busco/complete_busco_id.txt
sort busco/complete_busco_id.txt | uniq -c | awk '$1=="'$num_genomes'"{print $2}' > busco/final_busco_id.txt && rm busco/complete_busco_id.txt

cat busco/final_busco_id.txt | while read busco_id
do
	echo "Parsing..." $busco_id
	cat $FILE | while read assembly tol_id species genome
	do
		for faa in /lustre/scratch123/tol/projects/lepidoptera/data/insects/$species/analysis/$tol_id/busco/mm49_lepidoptera_odb10_metaeuk/run_lepidoptera_odb10/busco_sequences/single_copy_busco_sequences/$busco_id.faa
		do
			# Get the genomic coordinates of each single copy busco gene as it appears in the protein FASTA
			coordinates=`cat $faa | grep '>' | sed 's/:/\t/g' | sed 's/>//g' | sed 's/-/\t/g'`
			echo $coordinates $species $genome >> busco/$busco_id.tsv
		done
	done
	# We now need to change the contig/scaffold name because the one present in the protein FASTA file is based on the NCBI annotation
	python3 change_contig_name.py --species $FILE --busco busco/$busco_id.tsv && rm busco/$busco_id.tsv

	# We can now use the first species as a reference to extract from the cactus alignment an alignment in multiple alignment format (MAF)
	cat busco/$busco_id.bed | head -1 | while read contig start end length species genome;do hal2maf --refSequence $contig --start $start --length $length --refGenome $genome --onlyOrthologs --noAncestors $HAL busco/$busco_id.maf;done

	# Finally, we will check the consistency of the alignment by taking 100 bp upstream and downstream to have some flexibility
	reference=`cat busco/$busco_id.bed | head -1 | cut -f6`
	python3 check_consistency_single_copy_genes.py --maf busco/$busco_id.maf --busco busco/$busco_id.bed --bp 100 --refGenome $reference
done

