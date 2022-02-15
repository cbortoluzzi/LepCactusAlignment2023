#!/usr/bin/env python


# Author : @cb46


import argparse
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Assign number to contig/scaffold using NCBI assembly report')
parser.add_argument('--species', help = 'A tsv file with information on assembly, tol_id, species latin name, and name of the species as it appears in the cactus alignment')
parser.add_argument('--busco', help = 'A tab delimited file with genomic coordinates of single copy busco genes')



def parse_species(species_f):
	"""
	Read in the tsv file with a set of information on the species contained in the cactus alignment
	Input:
		species_f : tsv file
	Output:
		myspecies : dictonary with information on species latin name and species name as it appears in the alignment (key), tol_id, and assembly (value)
	"""
	myspecies = {}
	with open(species_f) as f:
		for line in f:
			assembly, tol_id, species_name, genome = line.strip().split()
			myspecies[species_name, genome] = [tol_id, assembly]
	return myspecies



def parse_assembly_report(myspecies):
	"""
	Read in the dictionary generated in the previous function and get the number of the contig/scaffold from the assembly report of each species, separately
	Input:
		myspecies : dictionary
	Output:
		mydict : a dcitionary with information on species latin name and species name as it appears in the alignment (key), and number of contig/scaffold (value)
	"""
	mydict = defaultdict(lambda: defaultdict(list))
	# Main path to all species
	path = '/lustre/scratch123/tol/projects/lepidoptera/data/insects/'
	for (name, genome), (tol_id, assembly) in myspecies.items():
		# Path to assembly report
		assembly_report = Path(path, name+'/assembly/release/'+tol_id+'/insdc/'+assembly+'_assembly_report.txt')
		with open(assembly_report) as f:
			for line in f:
				if not line.startswith('#'):
					line = line.strip().split()
					ncbi_contig = line[4]
					if line[1] == "assembled-molecule":
						chromosome_num = line[2]
						mydict[name, genome][ncbi_contig].append(chromosome_num.lstrip('0'))
					else:
						unplaced_scaffold = line[4]
						mydict[name, genome][ncbi_contig].append(unplaced_scaffold)
	return mydict



def change_ncbi_contig_to_chromosome_number(mydict, busco_gene, output):
	"""
	Change contig/scaffold name
	Input:
		mydict : dictionary
		busco_gene : tab delimited file with genomic coordinates of single copy busco genes
	Output:
		output : a tab delimited file with updated genomic coordinates of single copy busco genes
	"""
	with open(busco_gene) as f, open(output, 'a') as out:
		for line in f:
			contig, start, end, species, genome = line.strip().split()
			for (species, genome) in mydict:
				for ncbi_contig in mydict[species, genome][contig]:
					length = int(end) - int(start)
					out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(ncbi_contig, start, end, length, species, genome))




if __name__ == "__main__":
	args = parser.parse_args()
	species_information = parse_species(args.species)
	chromosome_number = parse_assembly_report(species_information)
	output = args.busco.replace('.tsv', '.bed')
	change_ncbi_contig = change_ncbi_contig_to_chromosome_number(chromosome_number, args.busco, output)

