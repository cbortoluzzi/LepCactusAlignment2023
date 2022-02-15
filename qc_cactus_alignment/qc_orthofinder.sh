#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
        echo "Usage: ./qc_orthofinder.sh <species.list> <cactus.alignment>"
        exit -1
fi


FILE=$1
HAL=$2



mkdir -p orthofinder && mkdir -p orthofinder/proteomes
cat $FILE | while read assembly tol_id species genome
do
        # Get the pep file of each species
        pep=`ls /lustre/scratch123/tol/projects/lepidoptera/data/insects/$species/analysis/$tol_id/gene/ensembl/latest/assembly.ensembl.pep.fa.gz`
        # Copy the pep file to the correct directory, while renaming it
        cp -r $pep orthofinder/proteomes/$species.pep.fa.gz
        # Decompress the pep file
        gunzip orthofinder/proteomes/$species.pep.fa.gz
        # Extract longest transcript per gene to speed up the analysis (we will alsoo try to run OrthoFinder on the total set of transcripts per gene)
        python /lustre/scratch123/tol/teams/durbin/users/cb46/softwares/OrthoFinder/tools/primary_transcript.py orthofinder/proteomes/$species.pep.fa
done

# Run OrthoFinder: we will first run it on the primary transcripts and then on the whole set of transcripts
time /lustre/scratch123/tol/teams/durbin/users/cb46/softwares/OrthoFinder/orthofinder -f orthofinder/proteomes/primary_transcripts/

# Get coordinates of single-copy orthogroups
python3 coordinates_orthogroups.py --single_copy orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Feb03/Orthogroups/Orthogroups_SingleCopyOrthologues.txt --orthogroups orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Feb03/Orthogroups/Orthogroups.tsv --pep orthofinder/proteomes/ --species species_list.tsv

# Check consistency of single copy orthogroups in the alignment
cat orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Feb03/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | while read orthogroup
do
        echo "Parsing ..." $orthogroup
        cat orthofinder/proteomes/orthogroups/$orthogroup.bed | head -1 | while read contig start end length species genome gene
        do
                hal2maf --refSequence $contig --start $start --length $length --refGenome $genome --onlyOrthologs --noAncestors $HAL orthofinder/proteomes/orthogroups/$orthogroup.maf
                python3 check_consistency.py --maf orthofinder/proteomes/orthogroups/$orthogroup.maf --bed orthofinder/proteomes/orthogroups/$orthogroup.bed --bp 100 --refGenome $genome
        done
done
