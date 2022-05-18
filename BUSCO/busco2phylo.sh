#!/bin/bash



# Author: @cb46



# Export path
export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/mafft-7.490-with-extensions/bin:$PATH
export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/trimal/source:$PATH
export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/standard-RAxML-master:$PATH



if [ -z $1 ]; then
	echo "Usage: ./busco2phylo.sh <species list file in tsv format> <name of busco input directory>"
	exit -1
fi



FILE=$1
INPUT_DIR=$2
DATE="$(date +'%Y%m%d')"
BUSCO="busco_"$DATE


mkdir -p $BUSCO


# Count number of genomes
num_genomes=`cat $FILE | wc -l`


# Obtain a list of Unique complete BUSCO genes, one per line
# The input file should look like this:
# GCA_904848185.1	fAcaLat1.1	fish	Acanthopagrus_latus
# GCA_900324465.3	fAnaTes1.3	fish	Anabas_testudineus
# where, in order, we have information on the assembly version, tol_id, class, and species_name
cat $FILE | while read assembly tol_id class species group
do
	cat $INPUT_DIR/$class/busco.v5/$species/vertebrata_odb10_metaeuk/run_vertebrata_odb10/full_table.tsv | grep -v '^#' | awk '$2=="Complete" {print $1}' >> $BUSCO/complete_busco_ids.txt
done
sort $BUSCO/complete_busco_ids.txt | uniq -c | awk '$1=="'$num_genomes'"{print $2}' > $BUSCO/final_busco_ids.txt && rm $BUSCO/complete_busco_ids.txt


# Obtain MAFFT alignment for each single-copy BUSCO gene
mkdir -p $BUSCO/mafft && mkdir -p $BUSCO/trimAl

cat $BUSCO/final_busco_ids.txt | while read busco_id
do
	mkdir -p $BUSCO/fasta/$busco_id
	cat $FILE | while read assembly tol_id class species group
	do
		for faa in $INPUT_DIR/$class/busco.v5/$species/vertebrata_odb10_metaeuk/run_vertebrata_odb10/busco_sequences/single_copy_busco_sequences/$busco_id.faa
		do
			# Reformat amino-acid fasta sequence of each single copy gene and species
			cat $faa | awk '/^>/{print ">'$species'"; next}{print}' > $BUSCO/fasta/$busco_id/$species\_$busco_id.fa
		done
	done

	cat $BUSCO/fasta/$busco_id/*.fa >> $BUSCO/mafft/$busco_id.aln

	# Perform alignment with MAFFT
	mafft --amino $BUSCO/mafft/$busco_id.aln > $BUSCO/mafft/$busco_id.aln.mafft

	Trim alignment with trimal
	trimal -in $BUSCO/mafft/$busco_id.aln.mafft -out $BUSCO/trimAl/$busco_id.aln.mafft.trimAl -automated1
done

# Generate supermatrix
python3 superalignment.py -i $BUSCO/trimAl -o $BUSCO/matrix

# Runn Raxml
raxmlHPC-PTHREADS-SSE3 -T 8 -f a -m PROTGAMMAJTT -N 100 -n my_busco_phylo -s busco_20220513/matrix/supermatrix.aln.faa -p 13432 -x 89090

