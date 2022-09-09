#!/usr/bin/env python



# Author : @cb46


import gzip
import argparse
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Generate multiple fasta file from a set of insect reference genome assemblies')
parser.add_argument('--i', help = 'Tab delimited file with information on species name, tol id, assembly id, and taxonomic order')
parser.add_argument('--p', help = 'Main path to reference genome assemblies')
parser.add_argument('--o', help = 'Output directory')


def get_reference_genome_assemblies(file, path, outDir):
	chrom_d = defaultdict(dict)
	fasta_d = {}
	with open(file) as f:
		for line in f:
			species_name, tol_id, assembly, order = line.strip().split()
			reference_genome = Path(path, species_name, 'assembly', 'release', tol_id, 'insdc', assembly + '.fasta.gz')
			assembly_report = Path(path, species_name, 'assembly', 'release', tol_id, 'insdc', assembly + '_assembly_report.txt')
			if reference_genome.is_file() and assembly_report.is_file():
				chromosomes = convert2chrom(assembly_report, chrom_d, species_name)
				fasta_d[species_name] = reference_genome

	output_f = Path(outDir, 'outgroups.fasta')
	with open(output_f, 'w') as output:
		for name, fasta in fasta_d.items():
			for record in SeqIO.parse(gzip.open(fasta, 'rt'), "fasta"):
				chrom = chromosomes[name][record.id]
				record.description = ''
				record.id = chrom + "_" + name
				SeqIO.write(record, output, 'fasta')


def convert2chrom(report_f, dictionary, species_name):
	with open(report_f) as f:
		for line in f:
			if not line.startswith('#'):
				lin = line.strip().split()
				if lin[1] == 'assembled-molecule':
					dictionary[species_name][lin[4]] = lin[2]
				elif lin[1] == 'unlocalized-scaffold':
					dictionary[species_name][lin[4]] = lin[0]
				else:
					dictionary[species_name][lin[4]] = lin[4]
	return dictionary




if __name__ == "__main__":
	args = parser.parse_args()
	p = Path(args.o)
	p.mkdir(parents=True, exist_ok=True)
	reference_genome_assemblies = get_reference_genome_assemblies(args.i, args.p, args.o)

  
