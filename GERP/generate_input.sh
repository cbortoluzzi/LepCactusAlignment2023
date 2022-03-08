#!/bin/bash



# Author: @cb46


export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/gffread:$PATH


if [ -z $1 ]; then
	echo "Usage: ./generate_input.sh <reference.genome> <cactus.alignment> <species.name> <tol.id>"
	exit -1
fi



REF=$1 # name of species to use as reference as it appears in the alignment (e.g. inachis_io_gca905147045v1)
HAL=$2 # name of cactus alignment in HAL format
SPECIES=$3 # species latin name (e.g. Inachis_io)
TOL=$4 # tol id of species (e.g. ilAglIoxx1.1)


mkdir -p sequences/$REF && mkdir -p annotation/$REF && mkdir -p assembly/$REF


# Print sequences of given genome in bed format
echo "Print sequences of given genome"
halStats --bedSequences $REF $HAL | cut -f1,3 > sequences/$REF/$REF.bed

# Extract, for each contig in the given genome, an alignment in MAF format
echo "Obtain MAF for each sequence"
cat sequences/$REF/$REF.bed | while read contig seq_length; do bsub -R'select[mem>12000] rusage[mem=12000]' -M12000 -q basement -n 15 -G rdgroup -J hal2maf -o output_%J -e error_%J hal2maf --refSequence $contig \
--refGenome $REF --noAncestors --onlyOrthologs $HAL sequences/$REF/$REF.$contig.maf; done

# Copy annotation (in gff3 format), assembly (in fasta format) and decompress
echo "Copy annotation and assembly"
cp /lustre/scratch123/tol/projects/lepidoptera/data/insects/$SPECIES/analysis/$TOL/gene/ensembl/latest/assembly.ensembl.genes.gff3.gz annotation/$REF/$REF.gff3.gz && gunzip annotation/$REF/$REF.gff3.gz
cp /lustre/scratch123/tol/projects/lepidoptera/data/insects/$SPECIES/assembly/release/$TOL/insdc/GCA*.fasta.gz assembly/$REF/$REF.fasta.gz && gunzip assembly/$REF/$REF.fasta.gz

# Change header of assembly. This needs to be done because contigs in the alignment are integers, while in the assembly are strings
echo "Change assembly nomenclature"
python3 change_fasta_header.py --fa assembly/$REF/$REF.fasta --report /lustre/scratch123/tol/projects/lepidoptera/data/insects/$SPECIES/assembly/release/$TOL/insdc/GCA*_assembly_report.txt --out assembly/$REF/$REF.renamed.fasta

# Remove transcripts with internal stop codons
echo "Run gffread to remove transcripts with internal stop codons"
gffread -v -V -g assembly/$REF/$REF.renamed.fasta annotation/$REF/$REF.gff3 -o annotation/$REF/$REF.filtered.gff3
# This shouldn't lead to anything

# Obtain merged CDS for all transcripts per gene
python3 merge_coding_sequence_per_transcript.py --gff3 annotation/$REF/$REF.gff3 --refGenome $REF --feature CDS

