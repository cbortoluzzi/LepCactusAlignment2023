#!/usr/bin/env python


# Author : @cb46


import os
import argparse
import subprocess
from ete3 import Tree
from Bio import AlignIO
from pathlib import Path
from collections import defaultdict



parser = argparse.ArgumentParser(description = 'Check consistency of single copy BUSCO genes genomic coordinates in cactus alignment')
parser.add_argument('--refGenome', help = 'Name of species to use as query')
parser.add_argument('--t', help = 'Phylogenetic tree obtained from the cactus alignment')
parser.add_argument('--f', help = 'Tab delimited species list file')
parser.add_argument('--b', help = 'List of complete, single copy BUSCO genes, one per line')
parser.add_argument('--hal', help = 'Input hal file')
parser.add_argument('--o', help = 'Output directory [default = busco_quality_check]', default = 'busco_quality_check')



main_path = '/lustre/scratch123/tol/teams/durbin/users/cb46/darwin'



def pairwise_comparisons(refGenome, tree, path):
	list_nodes = []
	tree = Tree(tree, format = 1)
	# Re-root the guide tree to the outgroup
	tree.set_outgroup('tinea_trinotella_gca905220615v1')
	for node in tree.traverse('postorder'):
		if node.is_leaf():
			# Remove the outgroup from the pairwise comparisons
			if node.name != refGenome:
				output_directory = Path(path, refGenome + '_vs_' + node.name)
				output_directory.mkdir(parents=True, exist_ok=True)
				list_nodes.append(node.name)
	mypair = list(map(lambda e: (e, refGenome), list_nodes))
	return mypair



def get_species_name_tol_id(species_list):
	species_d = {}
	with open(species_list) as f:
		next(f)
		for line in f:
			tol_id, pclass, order, superfamily, family, latin_name, assembly = line.strip().split()
			# Rename genome following the nomenclature in the cactus alignment
			genome = latin_name.lower() + '_' + assembly.replace('GCA_', 'gca').replace('.', 'v')
			species_d[genome] = [tol_id, assembly, latin_name]
	return species_d



def busco_genes(list_genes):
	busco_genes = []
	with open(list_genes) as f:
		for line in f:
			busco = line.strip()
			busco_genes.append(busco)
	return busco_genes



def change_assembly_coordinates(query, target, species_d):
	query_d, target_d = {}, {}
	# We need to change the chromosome nomenclature to match the one in the cactus alignment
	# For example, HG992306.1 is chromosome '1' in Tinea trinotella
	change_q = change_coordinates(main_path, query, species_d[query][0], species_d[query][1], species_d[query][2], query_d)
	change_t = change_coordinates(main_path, target, species_d[target][0], species_d[target][1], species_d[target][2], target_d)
	return change_q, change_t



def change_coordinates(main_path, genome, tol_id, assembly, species_name, dictionary):
	# Path to assembly report file: we will use this file to change the name of the chromosome
	assembly_report = Path(main_path, species_name, 'assembly/release', tol_id, 'insdc', assembly + '_assembly_report.txt')
	# Check that assembly report file exists
	if Path(assembly_report).is_file():
		with open(assembly_report) as f:
			for line in f:
				if not line.startswith('#'):
					line = line.strip().split()
					if line[1] == 'assembled-molecule':
						dictionary[line[4]] = line[2]
					else:
						dictionary[line[4]] = line[0]
	return dictionary



def get_busco_coordinates(query, target, species_d, busco_genes, change_q, change_t, hal, path):
	query_name = species_d[query][2]
	target_name = species_d[target][2]
	tol_id_query = species_d[query][0]
	tol_id_target = species_d[target][0]
	for gene in busco_genes:
		coordinates_q = Path(main_path, query_name, 'analysis', tol_id_query, 'busco/lepidoptera_odb10_metaeuk/run_lepidoptera_odb10/busco_sequences/single_copy_busco_sequences', gene + '.faa')
		coordinates_t = Path(main_path, target_name, 'analysis', tol_id_target, 'busco/lepidoptera_odb10_metaeuk/run_lepidoptera_odb10/busco_sequences/single_copy_busco_sequences', gene + '.faa')
		if coordinates_q.is_file() and coordinates_t.is_file():
			with open(coordinates_q) as gene_q, open(coordinates_t) as gene_t:
				for line_q, line_t in zip(gene_q, gene_t):
					if line_q.startswith('>') and line_t.startswith('>'):
						chrom_q = line_q.strip().split(':')[0].replace('>', '')
						chrom_t = line_t.strip().split(':')[0].replace('>', '')
						start_end_q = line_q.strip().split(':')[1]
						start_end_t = line_t.strip().split(':')[1]
						start_q, end_q = start_end_q.split('-')
						start_t, end_t = start_end_t.split('-')
						# Change chromosome same so that the one in the alignment corresponds to that from BUSCO
						chromosome_q = change_q[chrom_q]
						chromosome_t = change_t[chrom_t]
						# Obtain multiple sequence alignment for each gene of query and target species
						maf = multiple_alignment_format(query, chromosome_q, start_q, end_q, target, chromosome_t, start_t, end_t, gene, hal, path)




def multiple_alignment_format(query, chromosome_q, start_q, end_q, target, chromosome_t, start_t, end_t, gene, hal, path):
	mypos = defaultdict(list)
	len_q = int(end_q) - int(start_q)
	output_maf = Path(path, query + '_vs_' + target, gene + '.maf')
	output_qc = Path(path, query + '_vs_' + target, gene + '.qc')
	# Run hal2maf to obtain an alignment in multiple alignment format (MAF)
	cmd = 'hal2maf --refSequence %s --refGenome %s --start %s --length %s --noAncestors --onlyOrthologs --targetGenomes %s %s %s' %(chromosome_q, query, start_q, len_q, target, hal, output_maf)
	command = subprocess.check_output(cmd, shell = True).decode()
	# Parse multiple sequence alignment
	for multiple_alignment in AlignIO.parse(output_maf, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith(target):
				tgtChrom, tgtStart, tgtLen, tgtStrand, tgtSize = seqrec.id.split('.', 1)[1], seqrec.annotations['start'], seqrec.annotations["size"], seqrec.annotations["strand"], seqrec.annotations["srcSize"]
				# Change coordinates following the strand
				strand = ['+' if tgtStrand == 1  else '-']
				if strand[0] == '+':
					start = tgtStart
					end = tgtStart + tgtLen
				else:
					end = tgtSize - tgtStart
					start = tgtSize - (tgtStart + tgtLen)
				if tgtChrom == chromosome_t:
					mypos[tgtChrom, gene, start_t, end_t].append([start, end])
	with open(output_qc, 'w') as quality_check:
		for (chromosome, gene, start_g, end_g) in mypos:
			# Let's take 100 base pairs upstream and downstream the BUSCO coordinates to allow some flexiblity
			start_gene = int(start_g) - 100
			end_gene = int(end_g) + 100
			n = 0
			for (start, end) in mypos[chromosome, gene, start_g, end_g]:
				if start in range(start_gene, end_gene) and end in range(start_gene, end_gene):
					len = end - start
					n += len
			quality_check.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(query, chromosome_q, start_q, end_q, target, chromosome_t, start_g, end_g, n))
	# We can now remove the MAF file
	os.remove(output_maf)




if __name__ == "__main__":
	args = parser.parse_args()
	# Generate all possible pairwise comparisons between the query and the target species
	pairwise_combinations = pairwise_comparisons(args.refGenome, args.t, args.o)
	species_tol_id = get_species_name_tol_id(args.f)
	list_busco_genes = busco_genes(args.b)
	# Iterate over each pairwise comparison and each single-copy BUSCO gene
	for (target, query) in pairwise_combinations:
		# Check if folder is empty
		if len(os.listdir(directory) != 0:
			# Update genomic coordinates of BUSCO genes using the assembly report file
			assembly_coordinates_query, assembly_coordinates_target = change_assembly_coordinates(query, target, species_tol_id)
			# Obtain an alignment in multiple alignment format (MAF) for each BUSCO gene and check consistency
			evaluate_consistency = get_busco_coordinates(query, target, species_tol_id, list_busco_genes, assembly_coordinates_query, assembly_coordinates_target, args.hal, args.o)

