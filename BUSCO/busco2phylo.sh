#!/bin/bash



# Author: @cb46



# Export path
export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/mafft-7.490-with-extensions/bin:$PATH
export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/trimal/source:$PATH
export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/standard-RAxML-master:$PATH



if [ -z $1 ]; then
	echo "Usage: ./busco2phylo.sh <species list> <name of busco input directory>"
	exit -1
fi



file=$1
directory=$2


mkdir -p busco


# Count number of genomes
num_genomes=`cat $file | wc -l`
echo "We will build a phylogenetic tree on " $num_genomes "genomes"


# Obtain the complete list of single copy BUSCO genes shared by all species included in the input file
cat $file | while read assembly species_name
do
	if [ -e $directory/$species_name/vertebrata_odb10_metaeuk/run_vertebrata_odb10/full_table.tsv ]
	then
		printf '%s\t%s\n' "$directory/$species_name/vertebrata_odb10_metaeuk" "$species_name" >> path_to_busco.txt
	else
		printf '%s\t%s\n' "$directory/$species_name/$assembly/vertebrata_odb10_metaeuk" "$species_name" >> path_to_busco.txt
	fi
done

# Obtain unique single copy BUSCO genes
cat path_to_busco.txt | while read path species_name;do cat $path/run_vertebrata_odb10/full_table.tsv | grep -v '^#' | awk '$2=="Complete" {print $1}' >> busco/complete_busco_ids.txt;done
sort busco/complete_busco_ids.txt | uniq -c | awk '$1=="'$num_genomes'"{print $2}' > busco/final_busco_ids.txt && rm busco/complete_busco_ids.txt


# Obtain a MAFFT alignment for each single copy BUSCO gene
mkdir -p busco/mafft && mkdir -p busco/trimAl

cat busco/final_busco_ids.txt | while read busco_id
do
	mkdir -p busco/fasta/$busco_id
	cat path_to_busco.txt | while read path species_name
	do
		for faa in $path/run_vertebrata_odb10/busco_sequences/single_copy_busco_sequences/$busco_id.faa
		do
			# Reformat amino-acid fasta sequence of each single copy gene and species
			cat $faa | awk '/^>/{print ">'$species_name'"; next}{print}' > busco/fasta/$busco_id/$species_name\_$busco_id.fa
		done
	done

	cat busco/fasta/$busco_id/*.fa >> busco/mafft/$busco_id.aln
	# Perform alignment with MAFFT
	mafft --amino busco/mafft/$busco_id.aln > busco/mafft/$busco_id.aln.mafft
	rm busco/mafft/$busco_id.aln

	# Trim alignment with trimal
	trimal -in busco/mafft/$busco_id.aln.mafft -out busco/trimAl/$busco_id.aln.mafft.trimAl -automated1
done

# Generate supermatrix
python3 superalignment.py -i busco/trimAl -o busco/matrix

# Runn Raxml
raxmlHPC-PTHREADS-SSE3 -T 8 -f a -m PROTGAMMAJTT -N 100 -n my_busco_phylo -s busco/matrix/supermatrix.aln.faa -p 13432 -x 89090

