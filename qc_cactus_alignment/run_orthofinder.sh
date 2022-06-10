#!/bin/bash



# Author: @cb46


export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/OrthoFinder:$PATH


if [ -z $1 ]; then
        echo "Usage: ./run_orthofinder.sh <a tab delimited species list file>"
        exit -1
fi


FILE=$1


mkdir -p orthofinder/proteomes


cat $FILE | while read assembly tol_id class species_name superfamily
do
  # Copy proteome file 
  cp -r /lustre/scratch123/tol/projects/lepidoptera/data/insects/$species_name/analysis/$tol_id/gene/ensembl/latest/assembly.ensembl.pep.fa.gz orthofinder/proteomes/$species_name.pep.fa.gz
  # Unzip proteome file
  gunzip orthofinder/proteomes/$species_name.pep.fa.gz
  # Obtain primary transcript
  python3 primary_transcript.py --pep orthofinder/proteomes/$species_name.pep.fa
done

# Run OrthoFinder on all primary transcripts of all 49 lepidoptera species
orthofinder -f orthofinder/proteomes/primary_transcripts/

# Check consistency of single-copy orthogroups
python3 consistency_orthofinder.py --refGenome Tinea_trinotella --list_orthogroups orthofinder/proteomes/primary_transcripts/OrthoFinder/Results_Feb03/Orthogroups/Orthogroups_SingleCopyOrthologues.txt \
--species_list species_list.tsv --tree supermatrix_datafreeze_080621.treefile.pruned --hal Lepidoptera_88_way-202201.hal --o orthofinder_quality_check



