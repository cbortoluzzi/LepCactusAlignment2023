#!/bin/bash



# Author: @cb46


export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/OrthoFinder:$PATH


if [ -z $1 ]; then
        echo "Usage: ./run_orthofinder.sh"
        exit -1
fi


mkdir -p orthofinder/proteomes
cat species.tsv | while read assembly tol_id class species_name superfamily
do
  # Copy proteome file 
  cp -r /lustre/scratch123/tol/projects/lepidoptera/data/insects/$species_name/analysis/$tol_id/gene/ensembl/latest/assembly.ensembl.pep.fa.gz orthofinder/proteomes/$species_name.pep.fa.gz
  # Unzip proteome file
  gunzip orthofinder/proteomes/$species_name.pep.fa.gz
  # Obtain primary transcript
  python /lustre/scratch123/tol/teams/durbin/users/cb46/softwares/OrthoFinder/tools/primary_transcript.py orthofinder/proteomes/$species_name.pep.fa
done

# Run OrthoFinder on all primary transcripts of all 49 lepidoptera species
orthofinder -f orthofinder/proteomes/primary_transcripts/
