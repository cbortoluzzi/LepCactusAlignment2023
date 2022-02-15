#!/usr/bin/env python


# Author : @cb46


import argparse
from pathlib import Path
from collections import defaultdict




parser = argparse.ArgumentParser(description = 'Get genomic coordinates of single copy orthogroups (OrthoFinder)')
parser.add_argument('--single_copy', help = 'Single copy orthologues genes')
parser.add_argument('--orthogroups', help = 'Orthogroups')
parser.add_argument('--pep', help = 'Directory to all pep fasta files')
parser.add_argument('--species', help = 'A tab delimited file with information on assembly, tol_id, species latin name, and name as it appears in the cactus alignment')



def parse_species(species_f):
	mysp = {}
	with open(species_f) as f:
		for line in f:
			assembly, tol_id, species, genome = line.strip().split()
			mysp[species] = genome
	return mysp



def single_copy_orthologues_genes(single_copy):
	single_copy_dict = {}
	with open(single_copy) as f:
		for line in f:
			orthogroup = line.strip()
			single_copy_dict[orthogroup] = 'single-copy'
	return single_copy_dict



def orthogroup_per_species(orthogroups, single_copy_dict):
	myorthogroups = {}
	with open(orthogroups) as f:
		header = f.readline().split()[1:]
		for line in f:
			line = line.strip().split()
			orthogroup = line[0]
			if orthogroup in single_copy_dict.keys():
				species_gene = list(zip(line[1:], header))
				myorthogroups[orthogroup] =  species_gene
	return myorthogroups



def coordinates_pep_file(myorthogroups, files):
	mypep = {}
	mydict = defaultdict(list)
	for file in files:
		with open(file) as f:
			for line in f:
				if line.startswith('>'):
					header = line.strip().split()
					if len(header) == 7:
						gene = header[3].replace('gene:', '')
						tol_id, contig, start, end, strand = header[2].split(':')
						length = int(end) - int(start)
						mypep[gene] = [contig, start, end, length]
	for key, value in myorthogroups.items():
		for item in value:
			pep = item[0]
			coordinate = '\t'.join(map(str, mypep[pep]))
			mydict[key].append((item[1], pep, coordinate))
	return mydict



def save_to_file(mydict, path, mysp):
	for key in mydict:
		output_name = Path(path, key+'.bed')
		with open(output_name, 'w') as out:
			for item in mydict[key]:
				species = item[0].replace('.pep', '')
				out.write('{}\t{}\t{}\n'.format(item[2], species, mysp[species]))


				
if __name__ == "__main__":
	args = parser.parse_args()
	name_species_alignment = parse_species(args.species)
	single_copy = single_copy_orthologues_genes(args.single_copy)
	orthogroups = orthogroup_per_species(args.orthogroups, single_copy)
	list_pep = list(Path(args.pep).rglob('*.pep.fa'))
	genomic_coordinates = coordinates_pep_file(orthogroups, list_pep)
	path = Path(args.pep, 'orthogroups')
	path.mkdir(parents=True, exist_ok=True)
	to_file = save_to_file(genomic_coordinates, path, name_species_alignment)

	
