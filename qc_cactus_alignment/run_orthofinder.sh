#!/bin/bash



# Author: @cb46



export PATH=/lustre/scratch123/tol/teams/durbin/users/cb46/softwares/OrthoFinder:$PATH


mkdir -p orthofinder/proteomes


cat ../species_list.tsv | sed '1d' | while read tol_id class order superfamily family latin_name assembly
do
	proteome=/lustre/scratch123/tol/projects/lepidoptera/data/insects/$latin_name/analysis/$tol_id/gene/ensembl/latest/assembly.ensembl.pep.fa.gz
	if [[ -f "$proteome" ]]
	then
		# Copy proteome file
		cp -r /lustre/scratch123/tol/projects/lepidoptera/data/insects/$latin_name/analysis/$tol_id/gene/ensembl/latest/$assembly.*.pep.fa.gz orthofinder/proteomes/$latin_name.pep.fa.gz
		# Obtain primary transcript
		python3 primary_transcript.py --pep orthofinder/proteomes/$latin_name.pep.fa.gz --o orthofinder/proteomes/primary_transcripts
	fi
done

