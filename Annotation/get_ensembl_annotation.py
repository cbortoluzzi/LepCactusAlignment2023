#!/usr/bin/env python


# Author : @cb46


import os
import json
import argparse
from pathlib import Path


parser = argparse.ArgumentParser(description = 'Download assemblies from Ensembl Rapid Release')
parser.add_argument('-table', help = 'Ensembl table v.2021-11')
parser.add_argument('-metadata', help = 'Ensembl rapid release metadata')
parser.add_argument('-o', help = 'File with output directories, one per line')



def output_path(output):
	mypath = {}
	with open(output) as f:
		for line in f:
			taxon, path = line.strip().split()
			mypath[taxon] = path
	return mypath



def parse_table(table):
	mydict = {}
	with open(table) as f:
		next(f)
		for line in f:
			assembly_id, assembly_name, class_p, clade, order, family, species, size, contig_n50, scaffold_n50, assembly_level = line.strip().split('\t')
			species_name = species.replace(' ', '_')
			tol, assembly = assembly_id.split('.', 1)
			mydict[assembly] = [species_name, assembly_name, clade]
	return mydict



def download_genome_assembly(mydict, metadata, mypath):
	with open(metadata) as json_f:
		data = json.load(json_f)
		for i in data:
			assembly = str(i['assembly_accession'])
			geneset = str(i['geneset'])
			if assembly in mydict.keys():
				genome = mydict[assembly][0]+"-"+assembly+"-unmasked.fa.gz"
				ensembl_genome = '/'.join(["ftp://ftp.ensembl.org/pub/rapid-release/species", mydict[assembly][0], assembly, "genome", genome])
				ensembl_annotation = '/'.join(["ftp://ftp.ensembl.org/pub/rapid-release/species", mydict[assembly][0], assembly, "geneset", geneset, "*"])
				class_p = mydict[assembly][2]
				path = mypath[class_p]
				directory = Path(path, mydict[assembly][0])
				try:
					directory.mkdir(parents=True, exist_ok=False)
				except FileExistsError:
					print ("Folder is already there")
				cmd_genome = 'wget -nc %s --directory-prefix=%s' % (ensembl_genome, directory)
				cmd_annotation = 'wget -nc %s --directory-prefix=%s' % (ensembl_annotation, directory)
				print ("Downloading genome ....")
				os.system(cmd_genome)
				print ("Downloading annotation ...")
				os.system(cmd_annotation)



if __name__ == "__main__":
	args = parser.parse_args()
	output_f = output_path(args.o)
	parse_file = parse_table(args.table)
	download_assembly = download_genome_assembly(parse_file, args.metadata, output_f)

	
