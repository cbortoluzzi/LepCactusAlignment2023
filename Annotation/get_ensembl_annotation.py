#!/usr/bin/env python


# Author: @cb46


import os
import csv
import json
import argparse
from pathlib import Path



parser = argparse.ArgumentParser(description = 'Use Ensembl rapid release metadata to download annotations')
parser.add_argument('-table', help = 'Vertebrata table v.2021-11')
parser.add_argument('-csv', help = 'Ensembl rapid release species list')
parser.add_argument('-metadata', help = 'Ensembl rapid release metadata')



path = {'Fish' : '/lustre/scratch123/tol/teams/durbin/users/cb46/vertebrate_genomes_project/pan_vertebrates_alignment/fish/busco.v5',
	'Birds' : '/lustre/scratch123/tol/teams/durbin/users/cb46/vertebrate_genomes_project/pan_vertebrates_alignment/bird_alignment/busco.v5',
	'Mammals' : '/lustre/scratch123/tol/teams/durbin/users/cb46/vertebrate_genomes_project/pan_vertebrates_alignment/mammals/busco.v5',
	'Amphibians' : '/lustre/scratch123/tol/teams/durbin/users/cb46/vertebrate_genomes_project/pan_vertebrates_alignment/amphibians/busco.v5',
	'Reptiles' : '/lustre/scratch123/tol/teams/durbin/users/cb46/vertebrate_genomes_project/pan_vertebrates_alignment/reptiles/busco.v5',
	'Chordates' : '/lustre/scratch123/tol/teams/durbin/users/cb46/vertebrate_genomes_project/pan_vertebrates_alignment/chordates/busco.v5'
	}



def assert_assembly(table, ensembl):
	mydict = {}
	annotated = {}
	with open(table) as f:
		next(f)
		for line in f:
			line = line.strip().split('\t')
			assembly, assembly_name = line[0].split('.', 1)[1], line[1]
			mydict[assembly] = assembly_name


	with open(ensembl) as csv_f:
		csv_r = csv.reader(csv_f, delimiter = ',')
		next(csv_r)
		for line in csv_r:
			species_name, clade, assembly_name, assembly, annotation, rna_seq = line[0].replace(' ', '_'), line[3], line[5], line[6], line[9], line[11]
			if assembly in mydict.keys() and assembly_name == mydict[assembly]:
				annotated[assembly] = [species_name, assembly_name, clade, annotation, rna_seq]
	return annotated



def download_annotation(annotated, metadata, path):
	with open(metadata) as json_f:
		data = json.load(json_f)
		for i in data:
			assembly = str(i['assembly_accession'])
			if assembly in annotated.keys():
				geneset = str(i['geneset'])
				ensembl_rapid_release = '/'.join(['ftp://ftp.ensembl.org/pub/rapid-release/species', annotated[assembly][0], assembly, 'geneset', geneset, '*'])
				clade = annotated[assembly][2]
				directory = Path(path[clade], annotated[assembly][0])
				if Path(directory).is_dir():
					cmd = 'wget -nc %s --directory-prefix=%s' % (ensembl_rapid_release, directory)
					os.system(cmd)



if __name__ == "__main__":
	args = parser.parse_args()
	check_assembly = assert_assembly(args.table, args.csv)
	download = download_annotation(check_assembly, args.metadata, path)
  
  
