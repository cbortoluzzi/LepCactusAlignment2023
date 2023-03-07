#!/bin/bash



# Author: @cb46



export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/OrthoFinder:$PATH


mkdir -p orthofinder/proteomes


cat species_list.tsv | sed '1d' | while read tol_id class order superfamily family latin_name assembly
do
        proteome=/lustre/scratch123/tol/projects/lepidoptera/data/insects/$latin_name/analysis/$tol_id/gene/ensembl/latest/$assembly.ensembl.*.pep.fa.gz
        species_name=`echo $latin_name | tr '[:upper:]' '[:lower:]'`
        assembly_name=`echo $assembly | sed 's/\./v/g' | sed 's/_//g' | tr '[:upper:]' '[:lower:]'`
        genome=$species_name"_"$assembly_name
        # Copy proteome file
        cp -r /lustre/scratch123/tol/projects/lepidoptera/data/insects/$latin_name/analysis/$tol_id/gene/ensembl/latest/$assembly.*.pep.fa.gz orthofinder/proteomes/$genome.pep.fa.gz
        # Obtain primary transcript
        python3 primary_transcript.py --pep orthofinder/proteomes/$genome.pep.fa.gz --o orthofinder/proteomes/primary_transcripts
done

