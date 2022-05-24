#!/bin/bash



# Author: @cb46


export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/gffread:$PATH



if [ -z $1 ]; then
        echo "Usage: ./generate_input.sh <input hal file> <species name> <tol id> <name of reference genome>"
        exit -1
fi



HAL=$1
SPECIES=$2
TOL=$3
REF=$4



mkdir -p sequences/$REF && mkdir -p annotation/$REF && mkdir -p assembly/$REF


# Print sequences of given genome in bed format
halStats --bedSequences $REF $HAL | cut -f1,3 > sequences/$REF/$REF.bed


# Extract, for each contig in the given genome, an alignment in MAF format
cat sequences/$REF/$REF.bed | while read contig seq_length; do bsub -R'select[mem>18000] rusage[mem=18000]' -18000 -q basement -n 15 -G rdgroup -J hal2maf -o output_%J -e error_%J hal2maf --refSequence $contig --refGenome $REF \
--noAncestors --onlyOrthologs $HAL sequences/$REF/$REF.$contig.maf; done


# Copy annotation (in gff3 format), assembly (in fasta format) and decompress
cp /lustre/scratch123/tol/projects/lepidoptera/data/insects/$SPECIES/analysis/$TOL/gene/ensembl/latest/assembly.ensembl.genes.gff3.gz annotation/$REF/$REF.gff3.gz && gunzip annotation/$REF/$REF.gff3.gz
cp /lustre/scratch123/tol/projects/lepidoptera/data/insects/$SPECIES/assembly/release/$TOL/insdc/GCA*.fasta.gz assembly/$REF/$REF.fasta.gz && gunzip assembly/$REF/$REF.fasta.gz


# Change header of assembly. This needs to be done because contigs in the alignment are integers (e.g. 1, 2, 3, ...), while in the assembly are strings (e.g. HG996487.1, HG996488.1, HG996489.1, ...)
python3 change_fasta_header.py --fa assembly/$REF/$REF.fasta --report /lustre/scratch123/tol/projects/lepidoptera/data/insects/$SPECIES/assembly/release/$TOL/insdc/GCA*_assembly_report.txt --out assembly/$REF/$REF.renamed.fasta


# Remove transcripts with internal stop codons: this shouldn't lead to anything
gffread -v -V -g assembly/$REF/$REF.renamed.fasta annotation/$REF/$REF.gff3 -o annotation/$REF/$REF.filtered.gff3


# Obtain merged CDS of all transcripts per gene
python3 feature_selection.py --gff3 annotation/$REF/$REF.gff3 --refGenome $REF --feature CDS

